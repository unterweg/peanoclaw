PDT_PATH=../peano3/pdt

java -jar $PDT_PATH/lib/DaStGen.jar --inline --plugin PeanoHeapSnippetGenerator --naming PeanoHeapNameTranslator --include $PWD/../../src src/peanoclaw/dastgen/CellDescription.def src/peanoclaw/records/
java -jar $PDT_PATH/lib/DaStGen.jar --inline --plugin PeanoHeapSnippetGenerator --naming PeanoHeapNameTranslator --include $PWD/../../src src/peanoclaw/dastgen/VertexDescription.def src/peanoclaw/records/
java -jar $PDT_PATH/lib/DaStGen.jar --inline --plugin PeanoHeapSnippetGenerator --naming PeanoHeapNameTranslator --include $PWD/../../src src/peanoclaw/dastgen/Data.def src/peanoclaw/records/
java -jar $PDT_PATH/lib/DaStGen.jar --inline --plugin PeanoHeapSnippetGenerator --naming PeanoHeapNameTranslator --include $PWD/../../src src/peanoclaw/dastgen/PatchDescription.def src/peanoclaw/records/
java -jar $PDT_PATH/lib/DaStGen.jar --inline --plugin PeanoHeapSnippetGenerator --naming PeanoHeapNameTranslator --include $PWD/../../src src/peanoclaw/dastgen/LevelStatistics.def src/peanoclaw/statistics/
java -jar $PDT_PATH/lib/DaStGen.jar --inline --plugin PeanoHeapSnippetGenerator --naming PeanoHeapNameTranslator --include $PWD/../../src src/peanoclaw/dastgen/TimeIntervalStatistics.def src/peanoclaw/statistics/
