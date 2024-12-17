include common.mk

data/rmsd_*.csv: $(strip ${SCRIPT_DIR})/calc_rmsd.py
	python "$<" \
        --path-to-trajectory="$(strip ${TRAJECTORY_PATH})" \
        --path-to-trajectory-reference="$(strip ${TRAJECTORY_REFERENCE_PATH})" \
        --trajectory-stride="$(strip ${TRAJECTORY_STRIDE})" \
        --path-to-xray-reference="$(strip ${XRAY_REF})" \
        --filetype="$(strip ${FILETYPE})" \
        --pattern="$(strip ${PATTERN})" \
        --trajectory-length="$(strip ${TRAJECTORY_LENGTH})" \
        --dt-ns="$(strip ${DT_NS})" \
        --output-directory=data/