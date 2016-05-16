using FactCheck

include("detect_cycles.jl")

facts("Cycle Enumeration Test") do
    G = Graph()
    for i = 1:7
        addVertex(G, i)
    end
    addEdge(G, 1, 2)
    addEdge(G, 1, 3)
    addEdge(G, 2, 3)
    addEdge(G, 3, 4)
    addEdge(G, 4, 5)
    addEdge(G, 5, 6)
    addEdge(G, 6, 7)
    addEdge(G, 4, 7)

    cycles = OrderedSet{Int64}[]

    detectCycles(G, cycles)

    @fact length(cycles) --> 2

end
