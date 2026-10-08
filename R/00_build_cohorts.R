# Cohort construction from the four SEER*Stat case-listing exports (independent re-implementation of the authors'
# cohort script, extended columns). Reads NEC_DATA_DIR, writes patient-level cohort_*.csv and the aggregate
# eligibility flow flow_all.csv to NEC_WORK_DIR (never commit NEC_WORK_DIR: SEER Data Use Agreement).
# Run: Rscript R/00_build_cohorts.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))  # paths, NEC_BOOT, NEC_STRICT
suppressMessages(library(data.table))
DIR <- RAW_DIR
OUT <- DATA
COLS <- c("Histologic Type ICD-O-3","Behavior recode for analysis","Age recode with single ages and 90+","Primary Site - labeled",
  "Diagnostic Confirmation","Primary by international rules","Type of Reporting Source","Year of diagnosis","Survival Days","Survival months",
  "Record number recode","Patient ID","Combined Summary Stage with Expanded Regional Codes (2004+)","Summary stage 2000 (1998-2017)",
  "SEER cause-specific death classification","SEER other cause of death classification","Vital status recode (study cutoff used)","Sex",
  "Chemotherapy recode (yes, no/unk)","Radiation recode","RX Summ--Surg Prim Site (1998-2022)","RX Summ--Surg Prim Site 2023 (2023+)",
  "Survival months flag","Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)","Grade Recode (thru 2017)","Derived Summary Grade 2018 (2018+)",
  "Sequence number","Total number of in situ/malignant tumors for patient","Time from diagnosis to treatment in days recode",
  "SEER Combined Mets at DX-bone (2010+)","SEER Combined Mets at DX-brain (2010+)","SEER Combined Mets at DX-liver (2010+)","SEER Combined Mets at DX-lung (2010+)",
  "Mets at DX-Distant LN (2016+)","Mets at DX-Other (2016+)","Median household income inflation adj to 2024","Rural-Urban Continuum Code",
  "RX Summ--Systemic/Sur Seq (2007+)","Reason no cancer-directed surgery","Year of follow-up recode","COD to site recode")
build <- function(file, label, hist_keep, site_keep, diag_keep = "Positive histology") {
  r <- fread(paste0(DIR, file), sep="\t", colClasses="character", encoding="UTF-8", check.names=FALSE, na.strings=NULL, select=COLS)
  r[, row := .I]; flow <- list()
  st <- function(nm, x) { flow[[length(flow)+1]] <<- data.table(cohort=label, step=nm, records=nrow(x), patients=uniqueN(x[["Patient ID"]])) }
  st("exported", r)
  A <- r[trimws(`Histologic Type ICD-O-3`) %in% hist_keep]; st("histology", A)
  A <- A[`Behavior recode for analysis`=="Malignant"]; st("malignant", A)
  agenum <- suppressWarnings(as.numeric(sub("\\D.*$","",A[["Age recode with single ages and 90+"]])))
  A <- A[!is.na(agenum) & agenum>=18]; st("age>=18", A)
  A[, site5 := substr(`Primary Site - labeled`,1,5)]
  A <- A[site5 %in% site_keep]; st("site", A)
  A <- A[`Diagnostic Confirmation` %in% diag_keep]; st("diag_confirmation", A)
  A <- A[`Primary by international rules`=="Yes"]; st("intl_primary", A)
  A <- A[!(`Type of Reporting Source` %in% c("Autopsy only","Death certificate only"))]; st("not_autopsy_dco", A)
  A[, year := as.integer(`Year of diagnosis`)]
  A <- A[year>=2000 & year<=2023]; st("2000-2023", A)
  A[, sdays := suppressWarnings(as.numeric(`Survival Days`))]
  A <- A[!is.na(sdays) & sdays>=0]; st("survival_known", A)
  A[, recno := suppressWarnings(as.numeric(`Record number recode`))]
  setorder(A, `Patient ID`, year, recno, row)
  A <- A[!duplicated(`Patient ID`)]; st("one_per_patient", A)
  cs <- A[["Combined Summary Stage with Expanded Regional Codes (2004+)"]]; ss <- A[["Summary stage 2000 (1998-2017)"]]
  A[, stage := fifelse(year>=2004,
     fifelse(cs=="Localized only","Localized", fifelse(startsWith(cs,"Regional"),"Regional", fifelse(startsWith(cs,"Distant"),"Distant","Unknown"))),
     fifelse(ss=="Localized","Localized", fifelse(ss=="Regional","Regional", fifelse(ss=="Distant","Distant","Unknown"))))]
  A[, rectum := as.integer(site5 %in% c("C19.9","C20.9"))]
  A[, ev := fifelse(`SEER cause-specific death classification`=="Dead (attributable to this cancer dx)",1L,
             fifelse(`Vital status recode (study cutoff used)`=="Dead",2L,0L))]
  A[, unkcod := as.integer(`SEER cause-specific death classification`=="Dead (missing/unknown COD)")]
  A[, tm := pmax(sdays,0.5)/(365.25/12)]
  A[, age := as.numeric(sub("\\D.*$","",`Age recode with single ages and 90+`))]
  A[, male := as.integer(Sex=="Male")]
  A[, hist := trimws(`Histologic Type ICD-O-3`)]
  A[, chemo := as.integer(`Chemotherapy recode (yes, no/unk)`=="Yes")]
  A[, rad := as.integer(`Radiation recode` %in% c("Beam radiation","Radioactive implants (includes brachytherapy) (1988+)","Radioisotopes (1988+)","Combination of beam with implants or isotopes","Radiation, NOS  method or source not specified"))]
  A[, s1 := `RX Summ--Surg Prim Site (1998-2022)`][, s2 := `RX Summ--Surg Prim Site 2023 (2023+)`]
  A[, flagc := as.integer(`Survival months flag`=="Complete dates are available and there are more than 0 days of survival")]
  keep <- A[, .(pid=`Patient ID`, year, stage, site5, rectum, hist, ev, unkcod, tm, sdays, smonths=`Survival months`, flag=`Survival months flag`, age, male, chemo, rad, s1, s2, flagc,
     diag=`Diagnostic Confirmation`, race=`Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`, grade17=`Grade Recode (thru 2017)`, grade18=`Derived Summary Grade 2018 (2018+)`,
     seq=`Sequence number`, ntum=`Total number of in situ/malignant tumors for patient`, tt_days=`Time from diagnosis to treatment in days recode`,
     mets_bone=`SEER Combined Mets at DX-bone (2010+)`, mets_brain=`SEER Combined Mets at DX-brain (2010+)`, mets_liver=`SEER Combined Mets at DX-liver (2010+)`, mets_lung=`SEER Combined Mets at DX-lung (2010+)`,
     mets_dln=`Mets at DX-Distant LN (2016+)`, mets_other=`Mets at DX-Other (2016+)`, income=`Median household income inflation adj to 2024`, rural=`Rural-Urban Continuum Code`,
     sys_surg_seq=`RX Summ--Systemic/Sur Seq (2007+)`, no_surg_reason=`Reason no cancer-directed surgery`, fu_year=`Year of follow-up recode`, cod_site=`COD to site recode`,
     othercod=`SEER other cause of death classification`, vital=`Vital status recode (study cutoff used)`, reportsrc=`Type of Reporting Source`)]
  fwrite(keep, paste0(OUT,"cohort_",label,".csv"))
  cat(sprintf("%s: %d patients; stage: %s\n", label, nrow(keep), paste(names(table(keep$stage)), table(keep$stage), collapse=" ")))
  rbindlist(flow)
}
colo_sites <- c(sprintf("C18.%d", c(0,2:9)), "C19.9", "C20.9")
lung_sites <- sprintf("C34.%d", c(0,1,2,3,8,9))
panc_sites <- sprintf("C25.%d", c(0,1,2,3,4,7,8,9))
nen_all <- c("8013","8041","8042","8043","8044","8045","8240","8241","8242","8243","8244","8245","8246","8249")
ac_codes <- setdiff(c(as.character(8140:8389), as.character(8480:8490)), as.character(8240:8249))
sclc_all <- c("8002","8013","8040","8041","8042","8043","8044","8045")
panc_codes <- setdiff(c(as.character(8140:8389), as.character(8480:8490), "8500"), as.character(8240:8249))
cyt <- c("Positive histology","Positive exfoliative cytology, no positive histology")
f0 <- build(NEC_RAW_FILES[["NEN"]], "NEN_all", nen_all, colo_sites)        # all NEN codes; main NEC = 8013/8041-8045 subset
f1 <- build(NEC_RAW_FILES[["AC"]], "AC", ac_codes, colo_sites)
f2 <- build(NEC_RAW_FILES[["LUNG"]], "LUNG", sclc_all, lung_sites)
f3 <- build(NEC_RAW_FILES[["PANC"]], "PANC", panc_codes, panc_sites)
f4 <- build(NEC_RAW_FILES[["PANC"]], "PANC_cyt", panc_codes, panc_sites, cyt)
f5 <- build(NEC_RAW_FILES[["LUNG"]], "LUNG_cyt", sclc_all, lung_sites, cyt)
f6 <- build(NEC_RAW_FILES[["AC"]], "AC_cyt", ac_codes, colo_sites, cyt)
fl <- rbind(f0,f1,f2,f3,f4,f5,f6); fwrite(fl, paste0(OUT,"flow_all.csv")); print(fl)
