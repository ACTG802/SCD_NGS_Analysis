
#  iterate over all lines in the
#  requirements.txt file.

while read requirement; do conda install --yes $requirement; done < requirements_dev.txt