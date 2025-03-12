include common.mk

all: figure/rmsd.png

clean:
	rm -rf data
	rm -rf figure

figure/rmsd.png: $(strip ${SCRIPT_DIR})/plot_rmsd.py data/rmsd.csv
	python "$<" \
        --path-to-rmsd_dir=data/ \
        --output-directory=figure/

data/rmsd.csv: $(strip ${SCRIPT_DIR})/calc_rmsd.py
	python "$<" \
        --path-to-trajectory="$(strip ${TRAJECTORY_PATH})" \
        --path-to-trajectory-reference="$(strip ${TRAJECTORY_REFERENCE_PATH})" \
        --trajectory-stride="$(strip ${TRAJECTORY_STRIDE})" \
        --path-to-xray-reference="$(strip ${XRAY_REF})" \
        --trajectory-length="$(strip ${TRAJECTORY_LENGTH})" \
        --dt-ns="$(strip ${DT_NS})" \
        --output-directory=data/