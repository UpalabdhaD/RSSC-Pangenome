#!/usr/bin/env python3

import json
import sys


def merge(inputs):
    """ Merges the data in the given input handles and returns the result """
    assert len(inputs) > 1
    record_ids = set()
    merged = json.load(inputs[0])
    for record in merged["records"]:
        record_ids.add(record["id"])
    for other in inputs[1:]:
        data = json.load(other)
        assert data["schema"] == merged["schema"] == 1, "mismatching schema versions"
        # taxon wasn't enforced in schema 1, so use get instead
        assert data.get("taxon") == merged.get("taxon"), "mismatching taxons"
        for record in data["records"]:
            assert record["id"] not in record_ids, "duplicate record ids: %s" % record["id"]
            record_ids.add(record["id"])
            merged["records"].append(record)

    return merged


def _main(output_name, input_names):
    """ wrapper for catching exceptions and converting to exit values"""
    try:
        output = open(output_name, "w")
        inputs = [open(input_name) for input_name in input_names]
    except OSError as err:
        print(err, file=sys.stderr)
        return 1
    try:
        json.dump(merge(inputs), output)
    except json.decoder.JSONDecodeError:
        print("one or more files are not JSON formatted files", file=sys.stderr)
        return 1
    except KeyError as err:
        print("missing expected data key: %s" % err, file=sys.stderr)
        return 1
    except AssertionError as err:
        print(err, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: %s output input1 input2 [input3 input4 ...]" % sys.argv[0], file=sys.stderr)
        sys.exit(1)
    sys.exit(_main(sys.argv[1], sys.argv[2:]))
