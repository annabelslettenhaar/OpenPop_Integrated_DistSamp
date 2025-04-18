
library(tidyverse)
library(uuid)
library(lubridate)

## Convert old data format into 'new' without the unnecessary stuff

old_eve <- read.csv("data/ptar/event_data_1991_2014.csv")
old_occ <- read.csv("data/ptar/occurrence_data_1991_2014.csv")

## Event table
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


## Occurrence table
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
