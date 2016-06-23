#=======================================================
  AUTHOR: EMRE YAMANGIL, 2016
  Contact: emreyamangil@gmail.com, emreyamangil@lanl.gov
  Description:
  	Following code performs 
    1. a biconnected component decomposition of an undirected graph
    2. a DFS variant to find all cycles within each biconnected component. 
=======================================================#


using DataStructures

#=======================================================
  Class declarations
=======================================================#

type Vertex 
    label::Int64
    id::AbstractString
    AdjList::Set{Int64}
    visited::Bool
    included::Bool
    parent::Int64
    depth::Int64
    low::Int64
    pred::Int64
    function Vertex(x)
        v = new()
        v.label = x
        v.AdjList = Set{Int64}()
        v.visited = false
        v.included = false
        v.parent = -1
        return v
    end
end

type Edge
    from::Int64
    to::Int64
    function Edge(x, y)
        e = new()
        e.from = x
        e.to = y
        return e
    end
end

type Graph
    V::Vector{Vertex}
    E::Vector{Edge}
    function Graph()
        g = new()
        g.V = Vertex[]
        g.E = Edge[]
        return g
    end
end

function addEdge(G::Graph, x::Int64, y::Int64)
    push!(G.V[x].AdjList, y)
    push!(G.V[y].AdjList, x)
    push!(G.E, Edge(x,y))
end

function addVertex(G::Graph, x::Int64)
    push!(G.V, Vertex(x))
end

function min(a, b)
    return a <= b ? a : b
end

#=======================================================
  Biconnected component code
=======================================================#

function OutputComponent(S::Stack{Edge}, e::Edge, C::Vector{Vector{Edge}})
    Comp = Edge[]
    f = pop!(S)
    push!(Comp, f)
    while f.from != e.from || f.to != e.to
        f = pop!(S)
        push!(Comp, f)
    end
    push!(C, Comp)
    print("NEW COMPONENT: ")
    for i = 1:length(Comp)
        print("$(Comp[i]) ")
    end
    print("\n")
end

function DFS_visit(G::Graph, S::Stack{Edge}, u::Int64, depth::Int64, C::Vector{Vector{Edge}})
    depth += 1
    G.V[u].visited = true
    G.V[u].depth = depth
    G.V[u].low = depth
    for v in G.V[u].AdjList
        if !G.V[v].visited
            e = Edge(u,v)
            push!(S, e)
            G.V[v].parent = u
            DFS_visit(G, S, v, depth, C)
            if G.V[v].low >= G.V[u].depth
                OutputComponent(S, e, C)
            end
            G.V[u].low = min(G.V[u].low, G.V[v].low)
        elseif (G.V[u].parent != v) && (G.V[v].depth < G.V[u].depth)
            e = Edge(u,v)
            push!(S, e)
            G.V[u].low = min(G.V[u].low, G.V[v].depth)
        end
    end
end

function DFS(G::Graph, S::Stack{Edge}, source::Int64, C::Vector{Vector{Edge}})
    depth = 0
    N = length(G.V)
    for v = 1:N
        if !G.V[v].visited 
            DFS_visit(G, S, v, depth, C)
        end
    end

end

#=======================================================
  Depth first search code
=======================================================#

function GenerateSubGraph(C::Vector{Edge}, V::OrderedSet{Int64})
    for i = 1:length(C)
        push!(V, C[i].from)
        push!(V, C[i].to)
    end
end

function subGraphDFS_visit(G::Graph, E::Vector{Edge}, u::Int64, depth::Int64, V::OrderedSet{Int64})
    depth += 1
    G.V[u].visited = true
    G.V[u].depth = depth
    G.V[u].low = depth
    for v in G.V[u].AdjList
        if !G.V[v].visited && in(v, V)
            G.V[v].parent = u
            subGraphDFS_visit(G, E, v, depth, V)
        elseif G.V[u].parent != v && in(v, V)
            push!(E, Edge(u,v))
        end
    end
end

function subGraphDFS(G::Graph, E::Vector{Edge}, V::OrderedSet{Int64})
    depth = 0
    N = length(G.V)
    for i = 1:N
        G.V[i].visited = false
    end
    for i = 1:length(V)
        if !G.V[i].visited
            subGraphDFS_visit(G, E, i, depth, V)
        end
    end
end

function printGraph(G::Graph)
    for i = 1:length(G.V)
        println("Vertex $(G.V[i].label) has parent $(G.V[i].parent)")
    end
end

function printCycles(C::Vector{OrderedSet{Int64}})
    for i in 1:length(C)
        print("NEW CYCLE: (size $(length(C[i]))) ")
        for j in C[i].dict
            print("$j ")
        end 
        print("\n")
    end
end

function equalPath(p1::OrderedSet{Int64}, p2::OrderedSet{Int64})
    for x in p1
        if !in(x, p2)
            return false
        end
    end
    return true
end

function isNewPath(P::Vector{OrderedSet{Int64}}, path::OrderedSet{Int64})
    for i = 1:length(P)
        if equalPath(path, P[i])
            return false
        end
    end
    return true
end

function insertPath(P::Vector{OrderedSet{Int64}}, path_ini::OrderedSet{Int64})
    path = copy(path_ini)
    if isNewPath(P, path)
        push!(P, path)
      #  println("Found new path: $path")
    end
end

function findAllstPaths(G::Graph, P::Vector{OrderedSet{Int64}}, path_ini::OrderedSet{Int64}, u::Int64, t::Int64, V::OrderedSet{Int64})
    path = copy(path_ini)
    blah = G.V[u].AdjList
#    println("adjacent nodes to $u are $blah")
    
    for v in G.V[u].AdjList
        if v == t && length(path) == 1
            continue
        elseif !in(v, path) && in(v, V)
            if v == t
                push!(path, v)
                insertPath(P, path)
                delete!(path, v) 
            else
                push!(path, v)
               # println("current path: $path")
                findAllstPaths(G, P, path, v, t, V)
                delete!(path, v)
            end
        end
    end
end

#=======================================================
  Main method
=======================================================#

function detectCycles(G::Graph, cycles::Vector{OrderedSet{Int64}})

    C = Vector{Edge}[]
    S = Stack(Edge)
    source = 0
    DFS(G, S, source, C)

    println(length(C))
    for i = 1:length(C)
        println(i)
      
        V = OrderedSet{Int64}()
        GenerateSubGraph(C[i], V)
        E = Edge[]
        subGraphDFS(G, E, V)
       
        for i = 1:length(E)
            P = OrderedSet{Int64}()
            push!(P, E[i].from)
            findAllstPaths(G, cycles, P, E[i].from, E[i].to, V)
        end
 
    end

    printCycles(cycles)
end
