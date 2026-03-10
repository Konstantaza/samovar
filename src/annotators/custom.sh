custom_samovar () {
  local OPTIND opt
  while getopts "i:I:d:o:p:" opt; do
    case $opt in
      i) R1="$OPTARG" ;;
      I) R2="$OPTARG" ;;
      d) DB="$OPTARG" ;;
      o) OUT="$OPTARG" ;;
      p) PROPS="$OPTARG" ;;
    esac
  done

  concat=$(basename $R1 | sed 's/_.*//g')
  echo -e "\n" $concat "\n"
  mkdir -p $(dirname $OUT)
  echo -e "read_1\t562" > $OUT
  echo -e "read_2\t816" >> $OUT
  echo -e "read_3\t0" >> $OUT
}

custom_samovar "$@"