import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def matchEU(SE, SSIFULL):
    # Convert EventID to datetime format
    SSIFULL['EventID'] = pd.to_datetime(SSIFULL['EventID'].astype(str), format='%Y%m%d')

    # Initialize the new column in SE
    SE['EU'] = np.nan

    # Iterate over rows in SE
    for i, row in SE.iterrows():
        startTime = row['startTime']
        endTime = row['endTime']

        # Find rows in SSIFULL where EventID is within startTime and endTime
        rows = (SSIFULL['EventID'] >= startTime) & (SSIFULL['EventID'] <= endTime)

        # If any matching rows are found, take the sum of corresponding EU values
        if rows.any():
            SE.at[i, 'EU'] = SSIFULL.loc[rows, 'EU'].sum()
    
    return SE

def plotStormEnergyVsEU(SESSI):
    # Ensure there are no rows with NaN values
    SESSI = SESSI.dropna()

    # Normalize the stormEnergy and EU columns using min-max normalization
    SESSI['stormEnergy'] = (SESSI['stormEnergy'] - SESSI['stormEnergy'].min()) / (SESSI['stormEnergy'].max() - SESSI['stormEnergy'].min())
    SESSI['EU'] = (SESSI['EU'] - SESSI['EU'].min()) / (SESSI['EU'].max() - SESSI['EU'].min())

    # Create a scatter plot of normalized stormEnergy vs EU
    plt.figure()
    plt.scatter(SESSI['stormEnergy'], SESSI['EU'], c='blue', label='Data points')

    # Fit a linear regression model
    X = SESSI[['stormEnergy']]
    y = SESSI['EU']
    model = LinearRegression().fit(X, y)
    y_pred = model.predict(X)

    # Plot the regression line
    plt.plot(SESSI['stormEnergy'], y_pred, color='red', label='Linear fit')
    plt.title('Normalized Storm Energy vs EU')
    plt.xlabel('Normalized Storm Energy')
    plt.ylabel('Normalized EU')
    plt.grid(True)
    plt.legend()

    # Get the coefficients and R-squared value
    coeff = model.coef_[0]
    intercept = model.intercept_
    R2 = model.score(X, y)

    # Add the equation and the R-squared value to the plot
    equation = f'y = {coeff:.2f}x + {intercept:.2f}'
    Rsqrd = f'R² = {R2:.2f}'
    plt.text(0.05, 0.9, f'{equation}\n{Rsqrd}', transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

    plt.show()

def addNAOIndex(SESSI, naoDpM):
    # Initialize new column in SE
    SESSI['NAOIndex'] = np.nan

    # For each row in SE, find matching index_m value(s) in naoDpM
    for i, row in SESSI.iterrows():
        startTime = row['startTime']
        endTime = row['endTime']

        # Find rows in naoDpM where the date is within start and end times
        rows = (naoDpM['date'] >= startTime) & (naoDpM['date'] <= endTime)

        # If any matching rows were found, take the average index_m value
        if rows.any():
            SESSI.at[i, 'NAOIndex'] = naoDpM.loc[rows, 'index_m'].mean()
    
    return SESSI

# Example usage:
# SESSI = matchEU(SE, SSIFULL)
# plotStormEnergyVsEU(SESSI)
# SESSI = addNAOIndex(SESSI, naoDpM)

# For individual countries
# SE_AUT = addNAOIndex(SE_AUT, naoDpM);
# SE_CHE = addNAOIndex(SE_CHE, naoDpM);
# SE_CZE = addNAOIndex(SE_CZE, naoDpM);
# SE_DEU = addNAOIndex(SE_DEU, naoDpM);
# SE_DNK = addNAOIndex(SE_DNK, naoDpM);
# SE_EST = addNAOIndex(SE_EST, naoDpM);
# SE_FIN = addNAOIndex(SE_FIN, naoDpM);
# SE_FRA = addNAOIndex(SE_FRA, naoDpM);
# SE_GBR = addNAOIndex(SE_GBR, naoDpM);
# SE_HUN = addNAOIndex(SE_HUN, naoDpM);
# SE_IRL = addNAOIndex(SE_IRL, naoDpM);
# SE_LTU = addNAOIndex(SE_LTU, naoDpM);
# SE_LVA = addNAOIndex(SE_LVA, naoDpM);
# SE_NLD = addNAOIndex(SE_NLD, naoDpM);
# SE_NOR = addNAOIndex(SE_NOR, naoDpM);
# SE_POL = addNAOIndex(SE_POL, naoDpM);
# SE_SVK = addNAOIndex(SE_SVK, naoDpM);
# SE_SWE = addNAOIndex(SE_SWE, naoDpM);