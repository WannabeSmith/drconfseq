select
	sepsis3_fluids_w_subjectid.*,
	dod -- Date of death. Once combined with intime, compute survival
from
	(
		select
			subject_id,
			sepsis3_fluids_tbl.*
		from
			mimiciii.admissions as adm
			join (
				select
					total_fluids_tbl.hadm_id,
					fluids_24h_L,
					intime,
					age,
					gender,
					ethnicity,
					diabetes,
					bmi,
					first_service,
					elixhauser_hospital,
					sofa,
					lods,
					sirs,
					qsofa
				from
					public.sepsis3 as sepsis3 -- sepsis3 table created using alistairewj's scripts: https://github.com/alistairewj/sepsis3-mimic/blob/master/query/tbls/sepsis3.sql
					join (
						select
							iemv.hadm_id,
							sum(amount) / 1000 as fluids_24h_L -- Counting up the total fluids received in the first 24 hours of treatment
						from
							mimiciii.inputevents_mv as iemv -- Only considering patients in metavision following https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7036175/
							left join mimiciii.admissions as adm on iemv.hadm_id = adm.hadm_id
						where
							(
								-- The proxy we use for "is an inputevent a fluid?"
								ordercategoryname like '%Fluid%'
								or ordercategoryname like '%Crystalloids%'
								or ordercategoryname like '%Colloid%'
								or secondaryordercategoryname like '%Fluid%'
								or secondaryordercategoryname like '%Crystalloids%'
								or secondaryordercategoryname like '%Colloid%'
							)
							and amountuom = 'ml'
							and date_part(
								'day',
								iemv.starttime :: timestamp - adm.admittime :: timestamp
							) * 24 + date_part(
								'hour',
								iemv.starttime :: timestamp - adm.admittime :: timestamp
							) <= 24 -- Was the fluid administered within 24 hours of admission?
						group by
							iemv.hadm_id -- Count total fluids by admission id (this is also essentially by subject_id because we only consider the first hadm_id)
					) as total_fluids_tbl on sepsis3.hadm_id = total_fluids_tbl.hadm_id
				where
					excluded = 0 -- Apply exclusion criteria listed here https://github.com/alistairewj/sepsis3-mimic/blob/master/query/tbls/sepsis3.sql
					and suspected_infection_time_poe is not null
					and sofa >= 2 -- Using sepsis-3 definition according to Singer, Mervyn, et al. "The third international consensus definitions for sepsis and septic shock (Sepsis-3)." Jama 315.8 (2016): 801-810.
			) as sepsis3_fluids_tbl on sepsis3_fluids_tbl.hadm_id = adm.hadm_id
	) as sepsis3_fluids_w_subjectid -- table containing sepsis3 patients, 24hr fluid administered, with subject_id (for later joins)
	join mimiciii.patients as pt on pt.subject_id = sepsis3_fluids_w_subjectid.subject_id