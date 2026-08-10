# for metadata
datasets summary genome \
	taxon "Ralstonia solanacearum" \
	--assembly-level complete \
	--assembly-source genbank \
	--as-json-lines | \
	dataformat tsv genome > metadata_file.tsv

# for download (dryrun with --preview)
datasets download genome \
	taxon "Ralstonia solanacearum" \
	--assembly-source GenBank \
	--assembly-level complete \
	--include genome
datasets download genome taxon 'Ralstonia' --assembly-source GenBank --assembly-level complete --preview
