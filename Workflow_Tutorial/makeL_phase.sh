echo "date: $1"
echo "phase: $2" 
echo ${PWD}
date=$1
path=${PWD}
phase=$2
#ls ${path}/${date}/P*/LE/*net*.lc > ${date}/lc_LE.list
#ls ${path}/${date}/P*/ME/*net*.lc > ${date}/lc_ME.list
#ls ${path}/${date}/P*/HE/*25_80_net*.lc > ${date}/lc_HE.list

ls ${path}/${date}/P*/*${phase}_LE.pha > ${date}/LE_src_${phase}.list
ls ${path}/${date}/P*/1A*_LE_bkg.pha > ${date}/LE_bkg_${phase}.list
#ls ${path}/${date}/P*/selected_LE_bkg.pha > ${date}/LE_bkg_${phase}.list
ls ${path}/${date}/P*/*${phase}_LE.rsp > ${date}/LE_rsp_${phase}.list

ls ${path}/${date}/P*/*${phase}_ME.pha > ${date}/ME_src_${phase}.list
ls ${path}/${date}/P*/1A*_ME_bkg.pha > ${date}/ME_bkg_${phase}.list
#ls ${path}/${date}/P*/selected_ME_bkg.pha > ${date}/ME_bkg_${phase}.list
ls ${path}/${date}/P*/*${phase}_ME.rsp > ${date}/ME_rsp_${phase}.list

ls ${path}/${date}/P*/*${phase}_HE.pha > ${date}/HE_src_${phase}.list
ls ${path}/${date}/P*/1A*_HE_bkg.pha > ${date}/HE_bkg_${phase}.list
#ls ${path}/${date}/P*/selected_HE_bkg.pha > ${date}/HE_bkg_${phase}.list
ls ${path}/${date}/P*/*${phase}_HE.rsp > ${date}/HE_rsp_${phase}.list
