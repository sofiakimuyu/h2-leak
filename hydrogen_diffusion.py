def diffusion_time(volume_liters, diffusion_coefficient_cm2_per_s):
    # Convert volume from liters to cubic centimeters (cm^3)
    volume_cm3 = volume_liters * 1000

    # Calculate characteristic length
    length_cm = volume_cm3 ** (1 / 3)

    # Calculate time in seconds
    time_seconds = (length_cm ** 2) / (2 * diffusion_coefficient_cm2_per_s)

    return time_seconds


# Input values
volume_liters = 10  # total volume space in liters
diffusion_coefficient = 0.61  # diffusion coefficient of hydrogen in air at given conditions (cm^2/s)

time_sec = diffusion_time(volume_liters, diffusion_coefficient)

print(f"Estimated diffusion time: {time_sec:.1f} seconds (~{time_sec / 60:.2f} minutes)")
