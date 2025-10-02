import pandas as pd
import matplotlib.pyplot as plt


pathq = "/Users/matthijsnuus/Desktop/hymar/COMBA_Hymar/q_strain.tab"
path_aniso = "/Users/matthijsnuus/Desktop/hymar/COMBA_Hymar/q_strain_aniso.tab"
path_pp = "/Users/matthijsnuus/Desktop/hymar/COMBA_Hymar/pp_strain.tab"
pathw = "/Users/matthijsnuus/Desktop/hymar/winhausen_data/OPA-Z30-10MPa_OPA10.txt"


# read as whitespace-delimited, skip the first 2 metadata rows
df = pd.read_csv(
    pathq,
    sep=r"\s+",
    header=None,
    skiprows=2,
    engine="python"
)

# name columns (based on your filename/content)
df.columns = ["strain", "q"]
df['strain'] = df['strain'] * -1


# read as whitespace-delimited, skip the first 2 metadata rows
df_anis = pd.read_csv(
    path_aniso,
    sep=r"\s+",
    header=None,
    skiprows=2,
    engine="python"
)

# name columns (based on your filename/content)
df_anis.columns = ["strain", "q"]
df_anis['strain'] = df_anis['strain'] * -1

# read as whitespace-delimited, skip the first 2 metadata rows
df_pp = pd.read_csv(
    path_pp,
    sep=r"\s+",
    header=None,
    skiprows=2,
    engine="python"
)

# name columns (based on your filename/content)
df_pp.columns = ["strain", "pp"]
df_pp['strain'] = df_pp['strain'] * -1
df_pp['pp'] = df_pp['pp'] * 1e-6


df_z30_10 = pd.read_csv(
    pathw,
    sep="\t",          # it's tab-delimited
    header=0,          # first row is the header
    quotechar='"',     # headers are quoted
    encoding="utf-8",  # safe default
    engine="python"    # robust parser for mixed whitespace/quotes
)


# plot q vs strain
plt.figure()
plt.plot(df["strain"], df["q"], color='green', label='Modelled')
plt.plot(df_anis["strain"], df_anis["q"], color='red', label='Modelled-Anisotropy')

plt.plot(df_z30_10["axial strain[%]"], df_z30_10["diff. stress[MPa]"], color='blue', label='Measured')
plt.xlabel("Axial strain [%]")
plt.ylabel("q [MPa]")
plt.title("Z30 / 10 MPa")
plt.legend()
plt.tight_layout()
plt.show()


# plot q vs strain
plt.figure()
plt.plot(df_pp["strain"], df_pp["pp"] - 3.5 , color='green', label='Modelled')
plt.plot(df_z30_10["axial strain[%]"], df_z30_10["mean u[MPa]"] - 3.5 , color='blue', label='Measured')
plt.xlabel("Axial strain [%]")
plt.ylabel("excess pp [MPa]")
plt.title("Z30 / 10 MPa")
plt.legend()
plt.tight_layout()
plt.show()