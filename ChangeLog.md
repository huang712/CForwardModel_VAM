# 2020-09-20
- Add GPLv3 License

# 2020-09-11
- remove all functions and variables for Aaron's coordinates

# 2020-09-06
- add transformations for attitude dynamics (validated)
- use long int for quality_flags, double for utc_sec

The current CYGNSS antenna pattern (V6) is from version 2.1
We need the latest antenna pattern for version 3.0

# 2020-08-20
- remove checks on quality_flags in main.c

# 2020-08-05
- TxG = 1 in surface.c
- use Tx_eirp_watt instead of Tx_Power_dB in powerParm structure
- add sp_rx_gain in cygnss.c 
- correct the transformation from SPEC to Orbit (big change)
- index fixed in wind_interpolate (same as Thesis)

# 2020-07-22
- use CYGNSS reflection coefficient (cyg_R2)
- small update of CYGNSS GMF


