if (( $# < 2 )); then
    echo "Please provide two arguments:"
    echo "Folder to zip and the location of the writeup."
    exit 1
fi

rm "$1.zip" # deletes old zip file
zip -r "$1.zip" "$1" # zips the directory
zip -u "$1.zip" -j "$2" # adds the report to the zip file
