
import os
import pandas as pd
from aizynthfinder.reactiontree import ReactionTree
import argparse

class ImageRouteCreator:
    def __init__(self, fnJson, outputPath):
        self.fnJson = fnJson
        self.outputPath = outputPath


    def createImageRoutes(self):
        data = pd.read_json(self.fnJson, orient="table")
        data.to_json(self.fnJson, orient='table')
        print(f"Data saved to {self.fnJson}")

        if not data.empty and 'trees' in data.columns:
            all_trees = data['trees'].values
            if len(all_trees) > 0:
                trees_for_first_target = all_trees[0]
                if isinstance(trees_for_first_target, list):
                    for itree, tree in enumerate(trees_for_first_target):
                        if isinstance(tree, dict):
                            imagefile = os.path.join(self.outputPath, f"route{itree:03d}.png")
                            ReactionTree.from_dict(tree).to_image().save(imagefile)
                            print(f"Saved image to {imagefile}")
                        else:
                            print(f"Skipping non-dict tree at index {itree}")
                else:
                    print("Expected a list of trees, got:", type(trees_for_first_target))
            else:
                print("No trees found in the first target.")
        else:
            print("Data is empty or 'trees' column is missing.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create image routes from a JSON file.')
    parser.add_argument('inputJson', help='Path to the input JSON file.')
    parser.add_argument('outputDir', help='Path to the output directory.')

    args = parser.parse_args()

    creator = ImageRouteCreator(args.inputJson, args.outputDir)
    creator.createImageRoutes()

