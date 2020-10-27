# 2020-10-26
- Clean the code using Google style
- Minor changes to specular.c (use fabs in stead of abs)
- No longer calculate Katzberg model if GMF_OnOff==1 (wind.c and surface.c)

# 2020-09-20
- Add GPLv3 License

# 2020-09-11
- Remove all functions and variables for Aaron's coordinates

# 2020-09-06
- Add transformations for attitude dynamics (validated)
- Use long int for quality_flags, double for utc_sec

# 2020-08-20
- Remove checks on quality_flags in main.c

# 2020-08-05
- TxG = 1 in surface.c
- Use Tx_eirp_watt instead of Tx_Power_dB in powerParm structure
- Add sp_rx_gain in cygnss.c 
- Correct the transformation from SPEC to Orbit (big change)
- Index fixed in wind_interpolate (same as Thesis)

# 2020-07-22
- Use CYGNSS reflection coefficient (cyg_R2)
- Small update of CYGNSS GMF


