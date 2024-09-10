echo "date: $1"
echo "phase: $2"
echo ${PWD}
date=$1
phase=$2
cd ${PWD}/${date}
python ~/combine_HXMT_spec.py LE_src_${phase}.list LE_combine_${phase} -rmf_lst LE_rsp_${phase}.list -bkg_lst LE_bkg_${phase}.list
python ~/combine_HXMT_spec.py ME_src_${phase}.list ME_combine_${phase} -rmf_lst ME_rsp_${phase}.list -bkg_lst ME_bkg_${phase}.list
python ~/combine_HXMT_spec.py HE_src_${phase}.list HE_combine_${phase} -rmf_lst HE_rsp_${phase}.list -bkg_lst HE_bkg_${phase}.list
cd ${PWD}
