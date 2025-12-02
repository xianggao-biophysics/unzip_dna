import pandas as pd
import matplotlib.pyplot as plt

filename = 'unzipping_results.csv'
data = pd.read_csv(filename)
data = data.dropna()
extension_total = data['Extension_total(nm)']
force_average = data['Force_average(pN)']
plt.figure(figsize=(10, 6))
plt.plot(extension_total, force_average, label='Force vs Extension', color='blue')
plt.xlabel('Extension_total (nm)')
plt.ylabel('Force_average (pN)')
plt.title('Force vs Extension')
plt.legend()
plt.grid(True)
plt.show()
