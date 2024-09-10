echo "obsid: $1" 
echo "date: $2"
obsid=$1
date=$2
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0001 gti_tmp_${date}_0.0_0.1.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0102 gti_tmp_${date}_0.1_0.2.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0203 gti_tmp_${date}_0.2_0.3.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0304 gti_tmp_${date}_0.3_0.4.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0405 gti_tmp_${date}_0.4_0.5.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0506 gti_tmp_${date}_0.5_0.6.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0607 gti_tmp_${date}_0.6_0.7.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0708 gti_tmp_${date}_0.7_0.8.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0809 gti_tmp_${date}_0.8_0.9.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
python phase_resolved_HXMT_gti.py ${obsid}_phase ${obsid}_phase_0910 gti_tmp_${date}_0.9_1.0.txt -pipelineOutput /Volumes/WD_BLACK/Data/1A0535p262/HXMT/date_combine/${date}/${obsid}
