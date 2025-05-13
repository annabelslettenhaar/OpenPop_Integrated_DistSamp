
library(tidyverse)
library(uuid)
library(lubridate)
library(jsonlite)

## Convert old data format into 'new' without the unnecessary stuff

old_eve <- read.csv("data/ptar/event_data_1991_2014.csv")
old_occ <- read.csv("data/ptar/occurrence_data_1991_2014.csv")

#### Event table ####
new_eve <- old_eve %>%
  mutate(eventID = replicate(nrow(old_eve), UUIDgenerate()),
         transectlength = Lengde * 1000,
         stateProvince = case_when(
           OmradeNavn_New == "Åmotsdalen" ~ "Sør-Trøndelag",
           OmradeNavn_New == "Børgefjell" ~ "Nord-Trøndelag",
           OmradeNavn_New == "Møsvatn" ~ "Telemark",
           TRUE ~ NA_character_),
         municipality = case_when(
           OmradeNavn_New == "Åmotsdalen" ~ "Oppdal",
           OmradeNavn_New == "Børgefjell" ~ "Røyrvik",
           OmradeNavn_New == "Møsvatn" ~ "Tinn",
           TRUE ~ NA_character_),
         eventRemarks = "Line transect") %>%
  rename(area = OmradeNavn_New,
         locationID = LinNavn,
         year = Year) %>%
  select(eventID, area, locationID, year, transectlength, stateProvince, 
         municipality, eventRemarks)

# Now we need to add the 'Human observations' to the event table
# These describe the observation events, which are linked to the occurrence table with eventID
# Each eventID will have several occurrenceIDs if multiple types of birds are observed
# For example, it will provide the distance from the transect line that the birds were observed

lookup <- c("Åmotsdalen" = "Aamot", "Børgefjell" = "Bor", "Møsvatn" = "Mos")

# Remove observations of species we are not interested in
old_occ <- old_occ %>%
  filter(Art == 1)

new_occ <- old_occ %>%
  mutate(eventID = replicate(nrow(old_occ), UUIDgenerate()),
         eventRemarks = "Human observation",
         eventDate = as.Date(Dato, format = "%Y-%m-%d"),
         eventTime = format(ymd_hms(Kl, tz="UTC"), "%H:%M:%S"),
         stateProvince = case_when(
           Omrade == "Åmotsdalen" ~ "Sør-Trøndelag",
           Omrade == "Børgefjell" ~ "Nord-Trøndelag",
           Omrade == "Møsvatn" ~ "Telemark",
           TRUE ~ NA_character_),
         municipality = case_when(
           Omrade == "Åmotsdalen" ~ "Oppdal",
           Omrade == "Børgefjell" ~ "Røyrvik",
           Omrade == "Møsvatn" ~ "Tinn",
           TRUE ~ NA_character_),
         transectlength = NA_character_,
         locationID = paste0(lookup[Omrade], Linjenr_2)) %>%
  rename(area = Omrade,
         year = Aar,
         linedist = Linjeavst_2)

human_obs <- new_occ %>%
  select(eventID, area, locationID, eventDate, eventTime, year, stateProvince, 
         municipality, eventRemarks, linedist)

eve_complete <- bind_rows(human_obs, new_eve)

# Add parentEventID that contains the identifier of the corresponding transect line to human observations
eve_complete <- eve_complete %>%
  mutate(parentEventID = if_else(eventRemarks == "Line transect", eventID, NA_character_)) %>%
  group_by(locationID) %>%
  fill(parentEventID, .direction = "up") %>%
  ungroup() %>%
  mutate(parentEventID = if_else(eventRemarks == "Line transect", NA_character_, parentEventID))


#### Occurrence table ####
# Convert from wide to long format
new_occ_long <- new_occ %>%
  rename(Adult_Male = Adhann,
         Adult_Female = Adhunn,
         Adult_unknown = Adubest,
         Juvenile_unknown = Kyll_2,
         unknown_unknown = Ubest) %>%
  select(eventID, Adult_Male, Adult_Female, Adult_unknown,
         Juvenile_unknown, unknown_unknown) %>%
  pivot_longer(
    cols = -c(eventID),
    names_to = c("sex", "lifeStage"),
    names_sep = "_",
    values_to = "individualCount") %>%
  filter(individualCount > 0) %>%
  mutate(
    scientificName = "Lagopus lagopus",
    kingdom = "Animalia",
    phylum = "Chordata",
    class = "Aves",
    order = "Galliformes",
    genus = "Lagopus",
    specificEpithet = "lagopus",
    occurrenceID = replicate(nrow(new_occ_long), UUIDgenerate()))

write.csv(new_occ_long, file="data/ptar/occurrence_1991_2014.csv")  
write.csv(eve_complete, file="data/ptar/event_1991_2014.csv")  

#### Combine the old and new ptarmigan files into one ####

occ_old <- read.csv("data/ptar/occurrence_1991_2014.csv")
eve_old <- read.csv("data/ptar/event_1991_2014.csv")

occ_new <- read.csv("data/ptar/occurrence_2015_2020.csv")
eve_new <- read.csv("data/ptar/event_2015_2020.csv")

## Occurrence table

# First correct mistake in the old occ dataframe: swap sex and lifeStage column name
temp <- names(occ_old)[3]  
names(occ_old)[3] <- names(occ_old)[4]
names(occ_old)[4] <- temp

# Combine both dataframes
occ_total <- bind_rows(occ_new, occ_old)
occ_total <- occ_total %>%
  select(-X)

## Event table

eve_old <- eve_old %>%
  mutate(eventDate = ifelse(is.na(eventDate), as.Date("2001-01-01"), as.Date(eventDate))) %>%
  rename(verbatimLocality = area,
         sampleSizeValue = transectlength,
         oldDate = eventDate) # Dates look weird but will fix this later

# Change distance from transect line into the JSON format of the DwC layout
key <- "perpendicular distance in meters from transect line as reported by the field worker"

# Create a new column with the JSON strings in the old data
eve_old <- eve_old %>%
  mutate(dynamicProperties = sapply(linedist, function(x) {
    json_string <- paste0('{"', key, '":"', as.character(x), '"}')
    return(json_string)
    })
  ) %>%
  mutate(dynamicProperties = ifelse(grepl('"NA"', dynamicProperties), NA, dynamicProperties)) %>%
  select(-linedist, -X)

# Create a eventDate column with the correct year in the old data
# I do this because the 'year' column does not correspond to the eventDate
# After 2007, all eventDates are in 2019, so I assume 'year' is more accurate
eve_old$oldDate <- as.Date(eve_old$oldDate)
eve_old <- eve_old %>%
  mutate(eventDate = as.Date(paste(year, month(oldDate), day(oldDate), sep = "-"))) %>%
  select(-oldDate, - year)

eve_new$eventDate <- as.Date(eve_new$eventDate)
eve_total <- bind_rows(eve_new, eve_old)

eve_total <- eve_total %>%
  mutate(gyrArea = recode(verbatimLocality, 
                          "Kongsvoll" = "Dovrefjell",
                          "Møsvatn" = "Hardangervidda",
                          "Åmotsdalen" = "Dovrefjell",
                          "TOV-Åmotsdalen" = "Dovrefjell",
                          "TOV-Børgefjell" = "Børgefjell",
                          "TOV-Møsvatn" = "Hardangervidda"))

## Save the total files
write.csv(occ_total, file = "data/ptar/occurrence_total.csv")
write.csv(eve_total, file = "data/ptar/event_total.csv")


