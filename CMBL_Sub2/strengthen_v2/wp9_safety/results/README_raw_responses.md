# Saved openFDA responses

The responses that the scripts of this folder retrieved from the openFDA API are packed in `raw_responses.zip` (one JSON file per request; the file name is the SHA-1 of the request URL; each file holds the URL, the retrieval time and the response). They are public data. openFDA states that its data are unvalidated and must not be used for decisions about medical care.

To use the scripts, unpack the archive into a folder `raw/` next to it:

```
python -I -c "import zipfile; zipfile.ZipFile('raw_responses.zip').extractall('raw')"
```

`raw_responses_manifest.sha256` lists the SHA-256 of every file. The fingerprint recorded before the analysis (see `RAW_SHA256.txt` or the fingerprint lists) is the SHA-256 of the text made of one line `file name:SHA-256 of the file` for every file, sorted by file name; `freeze_raw.py` (work package 9) and `check_wp9_numbers.py` use this recipe.
