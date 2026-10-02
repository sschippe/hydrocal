    setDefaults("print summary: open", "Au75LMMDR.jac")
    grid       = Radial.Grid(Radial.Grid(false), rnt = 6.0e-6, h = 5.0e-2, hp = 6.0e-3, rbox = 10.0)
    drSettings = DielectronicRecombination.Settings(DielectronicRecombination.Settings(), multipoles = [E1,M1,E2], printBefore=true, 
                                                    gauges = [JAC.UseCoulomb, JAC.UseBabushkin])
    iConfigs   = [Configuration("1s^2 2s^2")] 
    eConfigs   = [Configuration("1s^2 2s 3s"), Configuration("1s^2 2s 3p"), Configuration("1s^2 2s 3d")] 
    fConfigs   = [Configuration("1s^2 2s^2 2p")] 
    addShells  = Basics.generateShellList(3, 3, [l for l=0:2]) 
    toShells   = Basics.generateShellList(2, 3, [0, 1, 2]) 
    finalConfigs  = Basics.generateConfigurations(fConfigs, [Shell("2s"), Shell("2p")] , toShells) 
    intermConfigs = Basics.generateConfigurationsWithAdditionalElectron(eConfigs, addShells) 
    intermConfigs = Basics.merge(intermConfigs,finalConfigs) 
    dcSettings = AsfSettings()   ## The standard version

    wa = Atomic.Computation(Atomic.Computation(), name="DR beryllium-like gold: 2s -> 3l", grid=grid, nuclearModel=Nuclear.Model(79.0),
                            initialConfigs      = [Configuration("1s^2 2s^2")],
                            intermediateConfigs = intermConfigs,
                            finalConfigs        = finalConfigs, finalAsfSettings = dcSettings,
                            processSettings     = drSettings )

    wb = perform(wa, output=true)
    setDefaults("print summary: close", "")
