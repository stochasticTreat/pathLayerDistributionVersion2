library('RUnit')

BiocGenerics:::testPackage("packageDir")

message("Running Tests")

test.suite <- defineTestSuite("example",
															dirs = system.file("unitTests", package = "packageDir"),
															testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)