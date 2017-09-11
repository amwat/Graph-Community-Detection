from igraph import *
from scipy import spatial


def qnewman(g, x, communityvertices):
    m = g.ecount()
    c = 1.0 / (2 * m)
    G = 0
    for i in communityvertices:
        if g.are_connected(i, x):
            G += g.es[g.get_eid(i, x)]["weight"]

    dx = sum([g.es[e]["weight"] for e in g.incident(x)])
    di = 0
    for v in communityvertices:
        di += sum([g.es[e]["weight"] for e in g.incident(v)])
    return c * (G - (dx * di * c))


def main():
    g = Graph()
    if len(sys.argv) != 2:
        print "usage : python sac1.py <alpha(0|0.5|1)>"
        exit(1)
    alpha = float(sys.argv[1])

    with open("data/fb_caltech_small_attrlist.csv") as f:
        attributes = next(f).strip().split(',')
        for line in f:
            attrs = map(float, line.strip().split(','))
            g.add_vertex(**dict(zip(attributes, attrs)))

    with open("data/fb_caltech_small_edgelist.txt") as f:
        for line in f:
            vertices = map(int, line.strip().split(' '))
            g.add_edge(vertices[0], vertices[1], **{"weight": 1.0})

    iteration = 0
    finalans = range(g.vcount())
    while iteration < 15:
        iteration += 1
        print "Iteration : ", iteration
        communities, community = phase1(g, alpha)
        phase2(g, community, communities, finalans)
        if len(community) == g.vcount():
            break
    finalcommunities = [[] for i in xrange(g.vcount())]
    for i in xrange(len(finalans)):
        finalcommunities[finalans[i]].append(i)
    print finalcommunities

    if alpha == 0:
        name = '0'
    elif alpha == 1:
        name = '1'
    else:
        name = '5'

    with open("communities_" + name + ".txt", "w") as f:
        for c in finalcommunities:
            f.write(','.join(map(str, c)))
            f.write('\n')


def phase2(g, community, communities, finalans):
    d = dict(enumerate(communities.keys()))
    d = {v: k for k, v in d.items()}
    finalans[:] = [d[community[x]] for x in finalans]
    community[:] = [d[x] for x in community]
    g.contract_vertices(community, combine_attrs=mean)
    g.simplify(combine_edges=sum)


def phase1(g, alpha):
    community = range(g.vcount())
    communities = {}
    for i in community:
        communities[i] = {i}
    cosinesim = [[0 for x in xrange(g.vcount())] for y in xrange(g.vcount())]
    for i in xrange(g.vcount()):
        for j in xrange(i, g.vcount()):
            cosinesim[i][j] = 1 - spatial.distance.cosine(g.vs[i].attributes().values(),
                                                          g.vs[j].attributes().values())
            cosinesim[j][i] = cosinesim[i][j]
    flag = True
    loop = 0
    while loop < 15 and flag:
        flag = False
        loop += 1
        print "\tLoop : ", loop
        for i in xrange(g.vcount()):
            maxgain = 0
            maxj = -1
            original = community[i]
            communities[original].discard(i)
            prevcos = []
            for x in communities[original]:
                prevcos.append(cosinesim[i][x])
            prevcosine = mean(prevcos)
            communities[original].add(i)
            for j in communities.keys():
                if original != j:
                    cossim = []
                    for x in communities[j]:
                        cossim.append(cosinesim[i][x])
                    community[i] = j
                    compositemodularitygain = alpha * qnewman(g, i, communities[j]) + (1 - alpha) * (mean(
                        cossim) - prevcosine)
                    if compositemodularitygain > maxgain:
                        maxgain = compositemodularitygain
                        maxj = j
            if maxj != -1:
                flag = True
                communities[original].discard(i)
                communities[maxj].add(i)
                community[i] = maxj
            else:
                community[i] = original

        for c in communities.keys():
            if not communities[c]:
                del communities[c]
        print "\t\tTotal : ", len(communities)
    return communities, community


if __name__ == "__main__":
    main()
