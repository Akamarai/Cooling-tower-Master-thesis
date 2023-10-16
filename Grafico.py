import matplotlib.pyplot as plt
import pandas as pd
import psychrolib

psychrolib.SetUnitSystem(psychrolib.SI)

df = pd.read_excel("meteo.xlsx")
# print(df.columns)

df.columns = ["time", "temperature", "humidity" , "pressure"]

df["time"] = pd.to_datetime(df["time"])
df["month"] = df["time"].dt.month

monthly_avg_temperature = df.groupby("month")["temperature"].mean()
monthly_avg_humidity = df.groupby("month")["humidity"].mean()
monthly_avg_pressure = df.groupby("month")["pressure"].mean()

# Calcola la temperatura a bulbo umido mensile
monthly_avg_WBT = []

for i in range(1, 13):
    DBT = monthly_avg_temperature[i]
    RH = monthly_avg_humidity[i] / 100
    P = monthly_avg_pressure[i] * 100

    WBT = psychrolib.GetTWetBulbFromRelHum(DBT, RH, P)
    monthly_avg_WBT.append(WBT)

month_names = [
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
]

plt.figure(figsize=(16, 9))

plt.plot(month_names, monthly_avg_temperature, label="Dry-bulb Temp [°C]", marker="o")
plt.plot(month_names, monthly_avg_humidity, label="Humidity ratio [%]", marker="o")
plt.plot(month_names, monthly_avg_WBT, label="Wet-bulb Temp [°C]", marker="o")

plt.xlabel("Month")
plt.ylabel("AVG value")
plt.title("Temperature and Humidity trends")
plt.legend()
plt.grid(True)

plt.show()
plt.savefig('andamento_medio_mensile.png')


data = {'Mese': month_names, 'Dry-bulb Temp [°C]': monthly_avg_temperature, 'Wet-bulb Temp [°C]': monthly_avg_WBT , 'Humidity ratio [%]': monthly_avg_humidity, 'Pressure [hPa]': monthly_avg_pressure}
table_df = pd.DataFrame(data)
table_df.to_excel('tabella_dati.xlsx', index=False)

