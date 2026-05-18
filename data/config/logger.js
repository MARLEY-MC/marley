{
  categories: {
    default: "info",
    physics.generator: "deBug",
    physics.generator.sampling: "TRACE",
  },

  out: [
    { stream: "stdout", max: "notice", min: "info" },
    { stream: "stderr", min: "warn" },
    { file: "full_log.txt", min: "trace" },
  ],
}
