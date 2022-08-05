rm(list=ls()); gc()
setwd("C:/repo/utils")

# install.packages("jsonlite")
# install.packages("stringdist")

library(jsonlite)
library(DBI)
library(magrittr)
library(tidyverse)
library(stringdist)

source("extract_util.R")

load_valueset.ncbo(vs_url="https://raw.githubusercontent.com/sxinger/PheCDM/main/valuesets/valueset_autogen/als-tx_output.json",
                   vs_name_str="riluzole")

load_valueset.curated(vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
                      vs_name_str = "myocardial infarction")

load_valueset.curated(vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
                      vs_name_str = "peripheral")

load_valueset.curated(vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
                      vs_name_str = "cerebrovascular")


load_valueset.curated(vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
                      vs_name_str = "bariatric surgery")

load_valueset.ecqm(vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_autogen/ecqm-allergy-intolerance.json",
                   vs_name_str = "statin")


load_valueset.curated(vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_autogen/ecqm-medication.json",
                      vs_name_str = "warfarin")


load_valueset(vs_template = "curated",
              vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
              vs_name_str = "myocardial infarction")


# load_valueset(vs_template = "curated",
#               vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
#               vs_name_str = "myocardial infarction",
#               dry_run = FALSE,
#               conn=myconn,
#               write_to_schema = "PUBLIC",
#               write_to_tbl = "TEMP")

source("extract_util.R")
collect_cdm(conn=NULL,
            patset = "PAT_ELIG",
            key_col_pat = c("PATID"),
            idx_dt_col = "INDEX_DATE",
            cdm_schema = "CMS_PCORNET_CDM",
            key_col_schema=c("PATID"),
            cdm_tbl = "DIAGNOSIS",
            ep_dt_col = "EVENT_DATE",
            vs_template = "curated",
            vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
            vs_name_str = "myocardial infarction",
            col_out = c("DX","DX_DATE"),
            write_back = FALSE,
            write_to_schema = "OBESITY",
            write_to_tbl = "TEMP",
            dry_run = TRUE)


source("extract_util.R")
load_valueset.curated(vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
                      vs_name_str = "bariatric surgery")
collect_cdm(conn=NULL,
            patset = "PAT_ELIG",
            key_col_pat = c("PATID"),
            idx_dt_col = "INDEX_DATE",
            cdm_schema = "CMS_PCORNET_CDM",
            key_col_schema=c("PATID"),
            cdm_tbl = "PROCEDURES",
            ep_dt_col = "EVENT_DATE",
            vs_template = "curated",
            vs_url = "https://raw.githubusercontent.com/sxinger/PheCDM/main/valusets/valueset_curated/vs-osa-comorb.json",
            vs_name_str = "bariatric surgery",
            col_out = c("PX","PX_DATE"),
            write_back = FALSE,
            write_to_schema = "OBESITY",
            write_to_tbl = "TEMP",
            dry_run = TRUE)

