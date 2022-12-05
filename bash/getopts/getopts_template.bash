#!/usr/bin/env bash


function usage() {
cat << EOF
usage: ${0} [OPTIONS] [ARGS]"

Version: ${VERSION}

OPTIONS:
    -h, --help    Show this message
    -d            Directory to scan [REQUIRED]
    --config      Path to a config file [REQUIRED]
    -v            display version

EOF
}

# single arg must be specified in 'optspec' and followed by ':' to accept values

optspec="hvd:-:"
while getopts "$optspec" optchar; do
    case "${optchar}" in
        -)
            case "${OPTARG}" in
                help)
                    usage
                    exit
                    ;;
                config)
                    CONFIG="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                *)
                    if [ "$OPTERR" = 1 ] && [ "${optspec:0:1}" != ":" ]; then
                        echo "Unknown option --${OPTARG}" >&2
                        usage
                        exit 1
                    fi
                    ;;
            esac;;
        d)
            SCAN_DIR=$OPTARG
            ;;
        v)
            version
            exit
            ;;
        h)
            usage
            exit
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

# Not mandatory but best practice to check argument
if [[ -z ${SCAN_DIR} ]] || [[ -z ${CONFIG} ]]
then
    log_error 'Options -d, -o, --config, -s are required'
    usage
    exit 1
fi

exit 0