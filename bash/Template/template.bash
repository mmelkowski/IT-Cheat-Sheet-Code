#!/usr/bin/env bash

VERSION='1.0.0'

function version() {
    echo ${VERSION}
}

function log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') [${0}] ${1}"
}


function log_info() {
    log "[INFO] ${1}"
}

function log_error() {
    log "[ERROR] ${1}" >&2
}

function usage() {
cat << EOF
usage: ${0} [OPTIONS] [ARGS]"

Version: ${VERSION}

OPTIONS:
    -h, --help    Show this message
    -d            Directory to scan (absolute path required) [REQUIRED]
    -o            Directory to output (absolute path required) [REQUIRED]
    --config      Path to a config file [REQUIRED]
    -v            display version

EOF
}

# single arg must be specified in 'optspec' and followed by ':' to accept values

optspec="hvd:o:-:"
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
        o)
            OUT_DIR=$OPTARG
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

if [[ -z ${SCAN_DIR} ]] || [[ -z ${OUT_DIR} ]] || [[ -z ${CONFIG} ]]
then
    log_error 'Options -d, -o, --config are required'
    usage
    exit 1
fi

# Load config file
source "${CONFIG}"

# Config Variables
folder="${SCAN_DIR}"
analysis_root_folder="${OUT_DIR}"
BIN_PATH="/usr/bin"

# Sub-bash example
nb_studies=$(find "${folder}" -maxdepth 1 -type f -name "*.tar.gz" | wc -l)

log_info "Scan ${nb_studies} studies in ${folder} ..."

# for loop
for study_archive_path in $(find "${folder}" -maxdepth 1 -type f -name "*.tar.gz" )
do
    study_archive=$(basename "${study_archive_path}")
    start_file="${study_archive_path}.ready"
    study_name=$(echo "${study_archive}" | sed 's/\(.*\)\.tar.gz/\1/')

    # if file exist
    if [[ -f "${start_file}" ]]
    then
        ${BIN_PATH}/prepare_analysis.py \
            --study_name "${study_name}" \
            --target_directory "${analysis_root_folder}"

        # Get command exit code
        exit_code=$?
        if [[ ${exit_code} -eq 0 ]]
        then
            log_info "Succes"
        else
            log_error "Fail"
            exit 1
        fi
    fi
done

exit 0