##########################################################################
# This file serves as the dataset registry for SciPy Datasets SubModule.
##########################################################################


data_repo_map = {
	"ascent": "dataset-ascent",
	"electrocardiogram": "dataset-ecg",
	"face": "dataset-face"
}

# To generate the SHA256 hash, use the command
# openssl sha256 <filename>
registry = {
	"ascent.dat": "e8a84939484463ab8051aedc5b40aa262ab33a91d6458a6cd13c6a1cad5a023d",
	"ecg.dat": "f20ad3365fb9b7f845d0e5c48b6fe67081377ee466c3a220b7f69f35c8958baf",
	"face.dat": "9d8b0b4d081313e2b485748c770472e5a95ed1738146883d84c7030493e82886"
}

registry_urls = {
    "ascent.dat": "https://raw.githubusercontent.com/scipy/dataset-ascent/main/ascent.dat",
    "ecg.dat": "https://raw.githubusercontent.com/scipy/dataset-ecg/main/ecg.dat",
    "face.dat": "https://raw.githubusercontent.com/scipy/dataset-face/main/face.dat"
}
