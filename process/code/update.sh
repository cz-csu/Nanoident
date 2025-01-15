
param=""

while getopts "p:" opt; do
  case $opt in
    p)
      param=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

destination="$param/home/nanodisco/code"

rsync -av --delete ./process/code "$destination"