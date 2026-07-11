{
  categories: {
    default: "trace",

    //// Physics — runtime event generation
    //physics.generator: "info",
    //physics.generator.sampling: "info",
    //physics.reaction: "info",
    //physics.reaction.xsec: "info",
    //physics.deexcitation: "info",
    //physics.deexcitation.hauser: "info",
    //physics.deexcitation.hauser.widths: "info",
    //physics.deexcitation.gamma: "info",
    //physics.coulomb: "info",
    physics.formfactor: "info",
    //physics.opticalmodel: "info",

    //// Initialization — configuration and data loading
    //init.config: "info",
    //init.config.source: "info",
    //init.config.target: "info",
    init.structure: "debug",
    //init.structure.decay: "info",
    //init.structure.masstable: "info",

    //// I/O — file management and output files
    io: "info",

    //// Application — CLI commands
    //app: "info",

    //// Tests
    //test: "info",
  },

  out: [
    { stream: "stdout", max: "notice", min: "info" },
    { stream: "stderr", min: "warn" },
    //{ file: "full_log.txt", min: "trace" },
  ],
}
