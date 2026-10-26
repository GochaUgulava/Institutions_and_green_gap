# ============================================================
# data.R  - data collection and master panel assembly
#  
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
  library(readxl)
  library(WDI)
})

# --------------------------- 1. Settings -------------------------------
PROJECT_DIR <- "G:/My Drive/Working Papers/Institutional Quality and Multidimensional Sustainable Development in EU Candidate Countries"
START_YEAR  <- 2009
END_YEAR    <- 2025
REF_YEAR    <- 2023
MAX_GAP     <- 2L
YEARS       <- START_YEAR:END_YEAR
EPI_YEAR    <- 2023L  

RAW_DIR       <- "data/raw"
PROCESSED_DIR <- "data/processed"
OUTPUT_FILE   <- file.path(PROCESSED_DIR, "master_panel_wide.csv")

WGI_FILE <- file.path(RAW_DIR, "wgidataset.xlsx")
VDEM_FILE <- file.path(RAW_DIR, "V-Dem-CY-Full+Others-v16.csv")
HDR_FILE <- file.path(RAW_DIR, "HDR25_Composite_indices_complete_time_series.csv")
EPI_FILE <- file.path(RAW_DIR, "epi2026.xlsx")
EPI_IND_DIR <- file.path(RAW_DIR, "epi_ind")
FH_FILE <- file.path(RAW_DIR, "freedom-score-fh.csv")
EIU_FILE <- file.path(RAW_DIR, "EIU_DI.csv")
OWID_CO2_FILE <- file.path(RAW_DIR, "owid-co2-data.csv")
OWID_ENERGY_FILE <- file.path(RAW_DIR, "owid_energy.csv")

if (dir.exists(PROJECT_DIR)) setwd(PROJECT_DIR)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------- 2. Countries and Indicators --------------------------- 
countries <- tribble(
  ~iso3, ~iso2, ~name, ~group,
  "AUT","AT","Austria","EU", "BEL","BE","Belgium","EU",
  "BGR","BG","Bulgaria","EU", "HRV","HR","Croatia","EU",
  "CYP","CY","Cyprus","EU", "CZE","CZ","Czechia","EU",
  "DNK","DK","Denmark","EU", "EST","EE","Estonia","EU",
  "FIN","FI","Finland","EU", "FRA","FR","France","EU",
  "DEU","DE","Germany","EU", "GRC","GR","Greece","EU",
  "HUN","HU","Hungary","EU", "IRL","IE","Ireland","EU",
  "ITA","IT","Italy","EU", "LVA","LV","Latvia","EU",
  "LTU","LT","Lithuania","EU", "LUX","LU","Luxembourg","EU",
  "MLT","MT","Malta","EU", "NLD","NL","Netherlands","EU",
  "POL","PL","Poland","EU", "PRT","PT","Portugal","EU",
  "ROU","RO","Romania","EU", "SVK","SK","Slovakia","EU",
  "SVN","SI","Slovenia","EU", "ESP","ES","Spain","EU",
  "SWE","SE","Sweden","EU",
  "ALB","AL","Albania","candidate",
  "BIH","BA","Bosnia and Herzegovina","candidate",
  "GEO","GE","Georgia","candidate", "MDA","MD","Moldova","candidate",
  "MKD","MK","North Macedonia","candidate", "MNE","ME","Montenegro","candidate",
  "SRB","RS","Serbia","candidate", "TUR","TR","Turkiye","candidate",
  "UKR","UA","Ukraine","candidate"
)

WDI_CODES <- c(
  gdp_pc = "NY.GDP.PCAP.PP.KD",
  gdp_ppp = "NY.GDP.MKTP.PP.KD",
  emp = "SL.EMP.TOTL.SP.ZS",
  life_exp = "SP.DYN.LE00.IN",
  pop = "SP.POP.TOTL",
  industry = "NV.IND.TOTL.ZS",
  pm25 = "EN.ATM.PM25.MC.M3"
)
WGI_CODES <- c("va", "pv", "ge", "rq", "rl", "cc")
EPI_CODES <- c(
  "PMD", "HFD", "NOD", "OZD", "UWD", "USD", "TBN", "TKP",
  "TCL", "RLI", "SNM", "PSU", "PBD", "PAR", "PAC", "KFL"
)
VDEM_CODES <- c(
  iso3           = "country_text_id",
  vdem_libdem    = "v2x_libdem",
  vdem_legcon    = "v2xlg_legcon",
  vdem_polyarchy = "v2x_polyarchy",
  vdem_frefair   = "v2xel_frefair",
  vdem_altinf    = "v2x_freexp_altinf",
  vdem_frassoc   = "v2x_frassoc_thick",
  vdem_jucon     = "v2x_jucon",
  vdem_cspart    = "v2x_cspart",
  vdem_corr      = "v2x_corr",
  vdem_rule      = "v2x_rule"
)

# --------------------------- 3. functions ------------------------------------
assert_unique_keys <- function(data, label) {
  if (anyDuplicated(data[c("iso3", "year")])) {
    stop("iso3-year key duplicates in the block: ", label, call. = FALSE)
  }
  data
}

fill_short_gaps <- function(x, max_gap = MAX_GAP) {
  if (sum(!is.na(x)) < 2L) return(x)
  x <- zoo::na.approx(x, maxgap = max_gap, na.rm = FALSE)
  zoo::na.locf(x, maxgap = max_gap, na.rm = FALSE)
}

read_wgi_sheet <- function(indicator) {
  data <- read_excel(WGI_FILE, sheet = indicator)
  code_col <- names(data)[str_detect(names(data), regex("^Economy \\(code\\)$", ignore_case = TRUE))]
  year_col <- names(data)[str_detect(names(data), regex("^Year$", ignore_case = TRUE))]
  estimate_col <- names(data)[str_detect(names(data), regex("^Governance estimate", ignore_case = TRUE))]

  if (length(code_col) != 1L || length(year_col) != 1L || length(estimate_col) != 1L) {
    stop(
      "Could not unambiguously identify the WGI columns on the sheet '", indicator,
      "'. Found: code=", paste(code_col, collapse = ","),
      "; year=", paste(year_col, collapse = ","),
      "; estimate=", paste(estimate_col, collapse = ","),
      call. = FALSE
    )
  }

  tibble(
    iso3 = as.character(data[[code_col]]),
    year = as.integer(data[[year_col]]),
    !!paste0("wgi_", indicator) := suppressWarnings(as.numeric(data[[estimate_col]]))
  ) %>%
    filter(iso3 %in% countries$iso3, year %in% YEARS)
}

read_epi_indicator <- function(code) {
  file_path <- file.path(EPI_IND_DIR, paste0(code, "_ind_na.csv"))
  if (!file.exists(file_path)) {
    warning("File EPI not find: ", file_path, call. = FALSE)
    return(NULL)
  }

  data <- read_csv(file_path, show_col_types = FALSE, na = c("NA", "")) %>%
    filter(iso %in% countries$iso3)
  value_cols <- setdiff(names(data), c("code", "iso", "country"))

  data %>%
    select(iso, all_of(value_cols)) %>%
    pivot_longer(-iso, names_to = "year_label", values_to = "value") %>%
    mutate(
      year = as.integer(str_extract(year_label, "\\d{4}")),
      code = code,
      value = if_else(value < 0, NA_real_, as.numeric(value))
    ) %>%
    select(iso3 = iso, year, code, value)
}

# --------------------------- 4. Panel --------------------------
panel_grid <- expand_grid(iso3 = countries$iso3, year = YEARS)

# WDI 
wdi_data <- WDI::WDI(
  country = countries$iso2,
  indicator = WDI_CODES,
  start = START_YEAR,
  end = END_YEAR,
  extra = FALSE
) %>%
  as_tibble() %>%
  left_join(select(countries, iso2, iso3), by = c("iso2c" = "iso2")) %>%
  select(iso3, year, all_of(names(WDI_CODES))) %>%
  filter(iso3 %in% countries$iso3) %>%
  assert_unique_keys("WDI")

# V-Dem
vdem_data <- read_csv(VDEM_FILE, show_col_types = FALSE, guess_max = Inf) %>%
  select(year, all_of(VDEM_CODES)) %>%
  mutate(year = as.integer(year)) %>%
  filter(iso3 %in% countries$iso3, year %in% YEARS) %>%
  assert_unique_keys("V-Dem")

# WID 
wid_data <- wid::download_wid(
  indicators = "sptinc",
  areas = countries$iso2,
  perc = "p0p50",
  ages = 992,
  pop = "j",
  years = YEARS
) %>%
  as_tibble() %>%
  left_join(select(countries, iso2, iso3), by = c("country" = "iso2")) %>%
  transmute(iso3, year = as.integer(year), p0p50 = as.numeric(value)) %>%
  filter(iso3 %in% countries$iso3, year %in% YEARS) %>%
  assert_unique_keys("WID")

# WGI
wgi_data <- map(WGI_CODES, read_wgi_sheet) %>%
  reduce(full_join, by = c("iso3", "year")) %>%
  assert_unique_keys("WGI")

# EPI 
epi_overall <- read_excel(EPI_FILE) %>%
  transmute(
    iso3 = str_to_upper(as.character(iso)),
    year = EPI_YEAR,
    epi = suppressWarnings(as.numeric(EPI.new))
  ) %>%
  filter(iso3 %in% countries$iso3) %>%
  distinct(iso3, year, .keep_all = TRUE) %>%
  assert_unique_keys("EPI overall")

epi_long <- map(EPI_CODES, read_epi_indicator) %>%
  keep(~ !is.null(.x)) %>%
  list_rbind()

wash_data <- epi_long %>%
  filter(code %in% c("UWD", "USD")) %>%
  pivot_wider(names_from = code, values_from = value) %>%
  mutate(
    value = case_when(
      !is.na(UWD) & !is.na(USD) ~ UWD * 0.6 + USD * 0.4,
      !is.na(UWD) ~ UWD,
      !is.na(USD) ~ USD,
      TRUE ~ NA_real_
    ),
    code = "WASH"
  ) %>%
  select(iso3, year, code, value)

epi_subindices <- epi_long %>%
  filter(!code %in% c("UWD", "USD")) %>%
  bind_rows(wash_data) %>%
  filter(year %in% YEARS) %>%
  mutate(variable = paste0("epi_", code)) %>%
  select(iso3, year, variable, value) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  assert_unique_keys("EPI subindices")

# UNDP HDR  
mys_data <- read_csv(HDR_FILE, show_col_types = FALSE) %>%
  filter(iso3 %in% countries$iso3) %>%
  select(iso3, matches("^mys_[0-9]{4}$")) %>%
  pivot_longer(-iso3, names_to = "year_label", names_prefix = "mys_", values_to = "mys") %>%
  transmute(iso3, year = as.integer(year_label), mys = as.numeric(mys)) %>%
  filter(year %in% YEARS) %>%
  assert_unique_keys("UNDP HDR")

# OWID energy and emissions
co2_data <- read_csv(OWID_CO2_FILE, show_col_types = FALSE) %>%
  transmute(
    iso3 = iso_code,
    year = as.integer(year),
    primary_energy = as.numeric(primary_energy_consumption),
    co2_mt = as.numeric(co2),
    co2_per_capita = as.numeric(co2_per_capita),
    consumption_co2 = as.numeric(consumption_co2),
    total_ghg = as.numeric(total_ghg)
  ) %>%
  filter(iso3 %in% countries$iso3, year %in% YEARS) %>%
  assert_unique_keys("OWID CO2")

energy_data <- read_csv(OWID_ENERGY_FILE, show_col_types = FALSE) %>%
  transmute(
    iso3 = iso_code,
    year = as.integer(year),
    renewables_share_elec = as.numeric(renewables_electricity_share_pct)
  ) %>%
  filter(iso3 %in% countries$iso3, year %in% YEARS) %>%
  assert_unique_keys("OWID energy")

# Freedom House and EIU indexes
fh_data <- read_csv(FH_FILE, show_col_types = FALSE) %>%
  transmute(
    iso3 = as.character(Code),
    year = as.integer(Year),
    fh_total = suppressWarnings(as.numeric(`Total democracy score`))
  ) %>%
  filter(iso3 %in% countries$iso3, year %in% YEARS) %>%
  assert_unique_keys("Freedom House")

eiu_data <- read_csv(EIU_FILE, show_col_types = FALSE) %>%
  transmute(
    iso3 = as.character(REF_AREA),
    year = as.integer(TIME_PERIOD),
    EIU_Index = as.numeric(OBS_VALUE),
    breakdown = COMP_BREAKDOWN_1
  ) %>%
  filter(
    iso3 %in% countries$iso3,
    year %in% YEARS,
    breakdown == "_T"
  ) %>%
  select(-breakdown) %>%
  assert_unique_keys("EIU Democracy Index")

# --------------------------- 5. Master panel ----------------------------
master_unfilled <- panel_grid %>%
  left_join(wdi_data, by = c("iso3", "year")) %>%
  left_join(vdem_data, by = c("iso3", "year")) %>%
  left_join(wid_data, by = c("iso3", "year")) %>%
  left_join(wgi_data, by = c("iso3", "year")) %>%
  left_join(epi_overall, by = c("iso3", "year")) %>%
  left_join(epi_subindices, by = c("iso3", "year")) %>%
  left_join(mys_data, by = c("iso3", "year")) %>%
  left_join(co2_data, by = c("iso3", "year")) %>%
  left_join(energy_data, by = c("iso3", "year")) %>%
  left_join(fh_data, by = c("iso3", "year")) %>%
  left_join(eiu_data, by = c("iso3", "year")) %>%
  left_join(select(countries, iso3, iso2, name, group), by = "iso3") %>%
  mutate(
    co2int = if_else(
      !is.na(co2_mt) & !is.na(gdp_ppp) & gdp_ppp > 0,
      co2_mt * 1e9 / gdp_ppp,
      NA_real_
    ),
    re_elec = renewables_share_elec
  ) %>%
  arrange(iso3, year)

write_csv(master_unfilled, file.path(PROCESSED_DIR, "master_panel_wide_unfilled.csv"), na = "")

master_panel <- master_unfilled %>%
  group_by(iso3) %>%
  mutate(across(
    where(is.numeric) & !all_of("year"),
    ~ fill_short_gaps(.x, MAX_GAP)
  )) %>%
  ungroup() %>%
  select(-any_of(c("re_share", "co2pc", "fossil_share_energy"))) %>%
  select(iso3, iso2, name, group, year, everything()) %>%
  arrange(iso3, year)


write_csv(master_panel, OUTPUT_FILE, na = "")

