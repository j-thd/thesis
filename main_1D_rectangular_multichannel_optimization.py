# File to configure optimization and run it

import optimization



if __name__ == "__main__":
    F_desired = 10e-3   # [N] Desired thrust
    P_max = 10          # [W] Maximum power consumption

    # Load all settings
    settings = optimization.settings.settings_1D_rectangular_multichannel
    optimization.run(
        F_desired=F_desired,
        P_max=P_max,
        settings=settings)