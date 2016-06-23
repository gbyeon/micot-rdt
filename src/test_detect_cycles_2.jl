include("detect_cycles.jl")

function create_random_graph(G,n,m)
	for j = 1:n
		addVertex(G, j)
	end
	for j = 1:m
		ind1 = convert(Int64, round((n-1)*rand())+1)
		ind2 = convert(Int64, round((n-1)*rand())+1)
		while ind2 in G.V[ind1].AdjList || ind1 == ind2
			ind1 = convert(Int64, round((n-1)*rand())+1)
			ind2 = convert(Int64, round((n-1)*rand())+1)
		end
		addEdge(G, ind1, ind2)
		addEdge(G, ind2, ind1)
	end
end

G = Graph()

create_random_graph(G, 10, 20)

cycles = OrderedSet{Int64}[]

detectCycles(G, cycles)
