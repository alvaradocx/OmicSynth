import streamlit as st

def app():
    st.title('Data Sources')

    st.text("Neuro disease constellation so far: \n \
    PD (Parkinson's) = Nalls 2019 \n \
    AD (Alzheimer's) = Jansen 2019 \n \
    ALS (Amyothrohic Lateral Sclerosis) = Nicolas 2018 \n \
    PSP (Progressive supranuclear palsy) = Hoglinger 2011 \n \
    FTD (Frontotemporal dementia) = Ferrari 2014")

    st.text("Intermediate so far: \n \
    MS (Multiple sclerosis) = IMSGC 2019 \n \
    T2D (Type 2 diabetes) = Mahajan 2018b \n \
    HTN (hypertension) = Neale lab 2021, from UKB mining \n \
    CVD (Cardiovascular disease) = Nelson 2017")

    st.text("Immuno disease constellation so far: \n \
    SLE (SLE) = Bentham 2015 \n \
    RA (Rheumatoid arthritis) = Stahl 2010 \n \
    CRDIBD (Crohn's disease + IBD) = Lange 2017 \n \
    CED (Celiac disease) = Dubois 2010")