'''
Early on in the covid-19 pandemic, it was recognized that the fatality rate in infected
individuals was much higher the older the patient was. The CDC released an estimate
of the fatality rate that you will use for this exercise:
Table 1: Fatalities per million infections:
0–17 years old: 20
18–49 years old: 500
50–64 years old: 6,000
65+ years old: 90,000
In order to determine the overall fatality rate, you must know the demographics of the
population. In this exercise, you will calculate the overall fatality rate for each country,
using the WorldDemographics.csv file provided. This file lists the number of individuals
that are living each each country by their age.
Your script should produced a csv file that has the following columns/rows.
Each row is a country
A column for the name of the country
A column for the total population of each country.
A column for the expected number of deaths if everyone in the country was infected by
covid-19.
A column for the percentage of the population that would die from covid-19.
To achieve this, I recommend you take the following steps:
1. Create a dataframe (dt1) using the Table 1 data where each row contains an age
(from 0-100) and expected fatality rate (e.g. 0.09).
2. Read in the WorldDemographics.csv file (dt2)
3. Join the dt1 and dt2 files together on Age.
4. Calculate the expected number of deaths for each age.
5. Group the data together by country to calculate the total population and expected
deaths for each country to create dt3.
6. Calculate a percentage died column in the dt3 table.
'''

import pandas as pd

# 1. Create a dataframe (dt1) using the Table 1 data where each row contains an age
# (from 0-100) and expected fatality rate (e.g. 0.09).
d = {'Age': [], 'FPM': []}
for i in range(101):
    d['Age'].append(i)
    if i < 18:
        d['FPM'].append(20/1000000)
    elif 18 <= i <= 49:
        d['FPM'].append(500/1000000)
    elif 50 <= i <= 64:
        d['FPM'].append(6000/1000000)
    else:
        d['FPM'].append(90000/1000000)

dt1 = pd.DataFrame.from_dict(d)

# 2. Read in the WorldDemographics.csv file (dt2)
dt2 = pd.read_csv('../hw/WorldDemographics.csv')
dt2 = dt2.drop(['Unnamed: 0'], axis=1)

# 3. Join the dt1 and dt2 files together on Age.
covid_dt = pd.merge(dt1, dt2, how='inner', on='Age')
# covid_dt = covid_dt.drop(['country_code', 'Level', 'continent_code'], axis=1)

# 4. Calculate the expected number of deaths for each age.
covid_dt['exp_death'] = covid_dt['FPM'] * covid_dt['#Alive']

# 5. Group the data together by country to calculate the total population and expected
# deaths for each country to create dt3.
dt3 = covid_dt.groupby('PopulationID').agg(total_polulation=('#Alive', 'sum'), total_exp_death=('exp_death', 'sum'))

# 6. Calculate a percentage died column in the dt3 table.
dt3['died_perc'] = dt3['total_exp_death']/dt3['total_polulation'] * 100

#print and save the dataframe as CSV
print(dt3)
dt3.to_csv('covid_fatality.csv', encoding='utf-8')