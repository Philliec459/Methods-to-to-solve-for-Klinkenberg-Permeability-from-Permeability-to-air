# Methods-to-to-solve-for-Klinkenberg-Permeability-from-Permeability-to-air

Methods to to solve for Klinkenberg Permeability from Permeability to air using Newton-Raphson and SciPy Optimization

Jones gives us in API Recommended Practices RP-40 that:

    b for Helium = 16.4*Kl**(-0.382) from equation A-7 of RP-40
    
    b for air   = 5.71*Kl**(-0.382) per RP-40 A-7 notes

By definition in equation 21 of RP-40:

    Kair = Kl (1 + b/p_mean)

    Kair = Kl + Kl*b_air/p_mean 

since:

    b_air = 5.71*Kl**(-0.382) derived from equation A-7 of RP-40, then

    Kair = Kl + Kl*5.71*Kl**(-0.382)/p_mean  

or 

    Kair =   Kl + 5.71*Kl**(1-0.382)/p_mean. 

Since,  

    0 = 5.71*Kl**0.618 + Kl * p_mean - Kair * p_mean ,
    
and the derivative of thsi equation is:

    f '(x) = 5.71 * 0.618 * Kl**(1 - 0.618) + p_mean

then we can iteratively solve for Kl knowing Kair and p_mean pressure. 


With this logic we can use an iterative approach to solve for Kl (Klinkenberg permeability) by rearranging the given equations and using an iterative method like Newton-Raphson to solve for Kl. 


```python
    # Given data
    K_air_mD = 46.6  # Measured permeability to air in mD
    p_mean_psi = 2.152  # Mean pressure in psi

    # Function to solve for Kl using Newton-Raphson method
    def calculate_Klinkenberg(K_air_mD, p_mean_psi):
        # Initial guess for Kl
        Kl = K_air_mD / 2  # Initial guess can be adjusted if needed

        # Define the function that needs to be solved to find Kl
        def f(Kl_guess):
            return 5.71 * (Kl_guess**(0.618)) + Kl_guess * p_mean_psi - K_air_mD * p_mean_psi

        # Define the derivative of the function for Newton-Raphson method
        def df(Kl_guess):
            return 5.71 * 0.618 * (Kl_guess**(-0.382)) + p_mean_psi

        # Set tolerance and maximum iterations for convergence
        tolerance = 1e-6
        max_iterations = 100
        
        iteration = 0
        while True:
            Kl_new = Kl - f(Kl) / df(Kl)

            if abs(Kl_new - Kl) < tolerance or iteration >= max_iterations:
                break

            Kl = Kl_new
            iteration += 1

        return Kl_new

    # Calculate Klinkenberg permeability using the function
    estimated_Kl = calculate_Klinkenberg(K_air_mD, p_mean_psi)
    print(f"The estimated Klinkenberg permeability is {estimated_Kl:.2f} mD.")
```

This code defines a function `calculate_Klinkenberg` that uses the Newton-Raphson method to iteratively solve for Kl given the provided equations and parameters. The initial guess for Kl is set to half the measured permeability to air, but this can be adjusted depending on the specific characteristics of the problem.

Running this code will output the estimated Klinkenberg permeability in milliDarcy (mD) based on the provided equations and iterative method. 

In the Newton-Raphson method, the derivative of the function with respect to the variable being solved for (in this case, Kl) is required to iteratively converge towards the solution. 

In the given equation:

    f(Kl) = 5.71 * Kl**{0.618} + Kl * p_mean - K_air * p_mean

The derivative of this function with respect to Kl can be calculated using standard differentiation rules. The derivative of the function f(Kl) with respect to Kl is denoted as f'(Kl) or df(Kl)/dKl. Applying the power rule and the derivative of a constant multiplied by a function, we get:

    f'(Kl) = df(Kl)/dKl = 5.71 * 0.618 * Kl**{0.618 - 1} + p_mean
    
    f'(Kl) = df(Kl)/dKl = 5.71 * 0.618 * Kl**{-0.382} + p_mean

This derivative, f'(Kl), is used in the Newton-Raphson method to update the guess for Kl iteratively. The derivative helps to find the rate of change of the function with respect to Kl at a particular point, which aids in determining how to update the guess to converge towards the root (or solution) of the equation.

In the provided Python code, the derivative function is defined as `df(Kl_guess)` and is used within the Newton-Raphson iteration to calculate the updated value of Kl. This iterative process continues until the difference between successive approximations falls below a specified tolerance or until a maximum number of iterations is reached.

The following uses the numbers supplied from Jones in RP-40 applying equation 31 where the Klinkenberg Permeability was 10.62. 
