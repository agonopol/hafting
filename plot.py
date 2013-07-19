from j2splat.visualize import *
import cjson

def main():
    frame = Frame("Grid Layer 2 Correlation")
    v = Visualize()
    with open("results.json") as f:
        json = cjson.decode(f.read())
        for l, trail in json.iteritems():
            resting = trail['resting']
            moving = trail['moving']
            for cell in resting.keys():
                v.append({"x":moving[cell], "y":resting[cell]})
    frame.append("Correlation Scatterplot", v.scatter())
    print frame.html()

if __name__ == "__main__":
    main()
