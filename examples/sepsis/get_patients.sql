select total_fluids_tbl.hadm_id, 
excluded, intime, dbsource, 
age, gender, ethnicity, diabetes, 
bmi, first_service, elixhauser_hospital, 
sofa, lods, sirs from public.sepsis3 as sepsis3
join
(select iemv.hadm_id, sum(amount) as total_fluids 
from mimiciii.inputevents_mv as iemv
left join mimiciii.admissions as adm
on iemv.hadm_id = adm.hadm_id
where (ordercategoryname like '%Fluid%' or
	   ordercategoryname like '%Crystalloids%' or
	   ordercategoryname like '%Colloid%' or
	   secondaryordercategoryname like '%Fluid%' or
	   secondaryordercategoryname like '%Crystalloids%' or
	   secondaryordercategoryname like '%Colloid%') and
	   amountuom = 'ml' and iemv.hadm_id = adm.hadm_id and
	   date_part('day', iemv.starttime::timestamp - adm.admittime::timestamp) * 24 + 
	   date_part('hour', iemv.starttime::timestamp - adm.admittime::timestamp) <= 24
	   group by iemv.hadm_id) as total_fluids_tbl
	   on sepsis3.hadm_id = total_fluids_tbl.hadm_id
	   where excluded = 0 and sepsis_cdc=1
