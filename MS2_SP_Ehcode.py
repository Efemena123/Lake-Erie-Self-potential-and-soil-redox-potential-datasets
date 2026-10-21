#%% my true inversion model

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import linregress
from tkinter import *
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
#import statsmodels.api as sm  # For linear regression


# Forward model block
def forward_model_poisson(charge_density, conductivity_profile, x_data, h):
    n = len(x_data)
    L = np.zeros((n, n))  # Solving kernel function with variable conductivity
    # Given constants
    #epsilon_0 = 8.85e-12  # Permittivity of free space (F/m)
    
    for i in range(n):
        if i > 0:
            L[i, i-1] = conductivity_profile[i] / (h**2)
        L[i, i] = -(conductivity_profile[i] + (conductivity_profile[i-1] if i > 0 else 0)) / (h**2)
        if i < n-1:
            L[i, i+1] = conductivity_profile[i] / (h**2) 
    
    predicted_sp = np.linalg.solve(L, - charge_density)
    return predicted_sp

# Inversion block
def inversion_function_poisson(charge_density, x_data, observed_sp, conductivity_profile, lambda_reg, h):
    predicted_sp = forward_model_poisson(charge_density, conductivity_profile, x_data, h)  # Pass conductivity_profile
    
    # Regularization to improve solution stability
    data_misfit = np.linalg.norm(observed_sp - predicted_sp)**2 #residual
    regularization = np.linalg.norm(np.diff(charge_density))**2
    return data_misfit + lambda_reg * regularization


# Define the core inversion and plotting functions
def run_inversion(h, lambda_reg=0.0001):
    try:
        # Use the globally loaded dataset
        global data
        x_data = np.array(data['Distance'])
        self_potential_data = np.array(data['SPc'])
        conductivity_profile = np.array(data['EC']) 
        Eh = np.array(data['Ehc'])
        SMC = np.array(data['SMC']) 
        T = np.array(data['T'])
        
        
        # Compute and store the residual
        '''data['SPC_SMC_residual'] = data['SPc'] - data['SMC']
        residual_sp = data['SPC_SMC_residual']
        print("Residual between SPc and SMC calculated and saved in dataframe.")
        '''
        # Linear Regression Model
        model = LinearRegression()
        model.fit(data[['Ehc']], data['SPc'])
        
        # Predict values
        y_predicted = model.predict(data[['Ehc']])
        
        # Calculate residuals
        residuals = data['SPc'] - y_predicted
        # Add residuals to DataFrame
        data['residuals'] = residuals
        
        
        n = len(x_data)

        # Initial guess and inversion
        initial_charge_density = np.ones(n)
        result = minimize(inversion_function_poisson, initial_charge_density,
                          args=(x_data, self_potential_data, conductivity_profile, lambda_reg, h),
                          method='L-BFGS-B', options={'maxiter': 100, 'disp': True})
        inverted_charge_density = result.x

        # Forward model with the inverted charge density
        predicted_sp = forward_model_poisson(inverted_charge_density, conductivity_profile, x_data, h)
        
        
       
        # Save results to CSV
        output_df = pd.DataFrame({
            'Distance': x_data,
            'Measured_SP': self_potential_data,
            'Inverted_SP': predicted_sp,
            'EC': conductivity_profile,
            'Inverted_Charge_Density': inverted_charge_density,
            'Soil_Eh': Eh,
            
            
        })

        output_df.to_csv("inverted_SP_results2.csv", index=False)
        print("Results saved as 'inverted_SP_results.csv'")


        # Call the function to display results
        display_results_in_gui(x_data, self_potential_data, predicted_sp, inverted_charge_density, Eh, SMC, T, residuals)
        
    except Exception as e:
        messagebox.showerror("Error", f"Failed to run inversion: {e}")
        return None

# Function to display the plots in the GUI
def display_results_in_gui(x_data, self_potential_data, predicted_sp, inverted_charge_density,Eh, SMC, T, residuals):
    fig, axs = plt.subplots(3, 2, figsize=(11, 9))

    # Plot Measured vs Inverted SP
    axs[0, 0].plot(x_data, self_potential_data, label='Measured SP')  # no line
    axs[0, 0].plot(x_data, predicted_sp , 'ro-', label='Inverted SP')         # no line
    axs[0, 0].set_title(f'Measure SP & Inverted SP (RMSE: {np.sqrt(np.mean((predicted_sp - self_potential_data)**2)):.2f} mV)')
    axs[0, 0].set_xlabel('Distance [m]')
    axs[0, 0].set_ylabel('Self-Potential (mV)')
    axs[0, 0].legend()

    # Plot inverted charge density
    axs[0, 1].plot(x_data, inverted_charge_density, 'ko-', label='Inverted current Density')
    axs[0, 1].set_title('Inverted current Density Distribution')
    axs[0, 1].set_xlabel('Distance [m]')
    axs[0, 1].set_ylabel('Current Density ($mA/m^2$)')
    axs[0, 1].legend()
    
    # Fit for Inverted SMC vs Eh
    slope_inverted, intercept_inverted, r_value_inverted, p_value_inverted, std_err_inverted = linregress(self_potential_data, predicted_sp)
    r_squared_inverted = r_value_inverted ** 2
        
    # Plot Eh vs SMC
    axs[1, 0].scatter(self_potential_data,predicted_sp, color='blue', label='Measured_SP vs Inversted_SP')
    axs[1, 0].plot(self_potential_data, intercept_inverted + slope_inverted * self_potential_data, 'k-', 
                   label=f'Fit: $R^2$={r_squared_inverted:.2f},P value={p_value_inverted:.1e},')
    axs[1, 0].set_xlabel('Measured_SP (mV)')
    axs[1, 0].set_ylabel('Inverted_SP (mV)')
    axs[1, 0].set_title('Measured_SP vs Inversted_SP')
    axs[1, 0].legend()

    # Fit for Inverted SP vs Eh
    slope_inverted, intercept_inverted, r_value_inverted, p_value_inverted, std_err_inverted = linregress(Eh, predicted_sp-1.83)
    r_squared_inverted = r_value_inverted ** 2

    # Plot Inverted SP vs Eh
    axs[1, 1].scatter(Eh, predicted_sp-1.83, color='black', label='Residual SP vs Eh')
    axs[1, 1].plot(Eh, intercept_inverted + slope_inverted * Eh, 'k-', label=f'Fit: $R^2$={r_squared_inverted:.2f}, P value={p_value_inverted:.1e}')
    axs[1, 1].set_xlabel('Soil Eh [mV]')
    axs[1, 1].set_ylabel('Residual SP (mV)')
    axs[1, 1].set_title('Residual SP vs Eh')
    axs[1, 1].legend()
    
      
 
    slope_inverted, intercept_inverted, r_value_inverted, p_value_inverted, std_err_inverted = linregress(SMC, Eh)
    r_squared_inverted = r_value_inverted ** 2
        
    # Plot Eh vs SMC
    axs[2, 0].scatter(SMC, Eh,color='blue', label='Eh vs SMC')
    axs[2, 0].plot(SMC, intercept_inverted + slope_inverted * SMC, 'k-', label=f'Fit: $R^2$={r_squared_inverted:.2f},P value={p_value_inverted:.2e}')
    axs[2, 0].set_xlabel('SMC (%)')
    axs[2, 0].set_ylabel('Soil Eh (mV)')
    axs[2, 0].set_title('Soil Eh vs SMC')
    axs[2, 0].legend()
    
    
 
    slope_inverted, intercept_inverted, r_value_inverted, p_value_inverted, std_err_inverted = linregress(Eh, self_potential_data)
    r_squared_inverted = r_value_inverted ** 2
        
    # Plot Eh vs SMC
    axs[2, 1].scatter(Eh, self_potential_data,color='blue', label='Eh vs measure SP')
    axs[2, 1].plot(Eh, intercept_inverted + slope_inverted * Eh, 'k-', label=f'Fit: $R^2$={r_squared_inverted:.2f},P value={p_value_inverted:.2e}')
    axs[2, 1].set_xlabel('Eh (mV')
    axs[2, 1].set_ylabel('measured SP (mV)')
    axs[2, 1].set_title('Soil Eh vs self_potential_data')
    axs[2, 1].legend()
    
    fig.tight_layout()
    plt.show()
    
# Function to load the file and initialize the inversion parameters
def load_file():
    global data
    file_path = filedialog.askopenfilename()
    if file_path:
        data = pd.read_csv(file_path)
        messagebox.showinfo("Success", "Dataset loaded successfully!")
    else:
        messagebox.showerror("Error", "Failed to load dataset.")

# Function to enter the value for h
def enter_h_value():
    global h
    try:
        h = float(h_entry.get())
        messagebox.showinfo("Success", f"h value set to {h}")
    except ValueError:
        messagebox.showerror("Error", "Please enter a valid numeric value for h.")

# Setup the Tkinter GUI
root = Tk()
root.title("SP Inversion Tool")
root.geometry("800x600")

# Create frames for buttons and plots
button_frame = Frame(root)
button_frame.pack(side=TOP, fill=X)

plot_frame = Frame(root)
plot_frame.pack(side=TOP, fill=BOTH, expand=True)

# Create buttons for forward model, inversion, and exit
label = Label(button_frame, text="SP Data Inversion", font=("Arial", 16))
label.pack(pady=10)

# Load file button
load_button = Button(button_frame, text="Load Dataset", command=load_file)
load_button.pack(pady=5)

# Entry for h value
h_label = Label(button_frame, text="Enter h value (step size for discretization):")
h_label.pack(pady=5)
h_entry = Entry(button_frame)
h_entry.pack(pady=5)
h_button = Button(button_frame, text="Set h Value", command=enter_h_value)
h_button.pack(pady=5)

# Run inversion button
run_button = Button(button_frame, text="Run Inversion", command=lambda: run_inversion(h))
run_button.pack(pady=10)

# Run the GUI loop
root.mainloop()



