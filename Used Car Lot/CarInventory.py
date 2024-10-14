from CarInventoryNode import CarInventoryNode
from Car import Car

class CarInventory:

    def __init__(self):
        self.root = None

    def addCar(self, car):

        if self.root:
            self._put(car, self.root)
        else:
            self.root = CarInventoryNode(car)

    def _put(self, car, currentNode):
        """Helper function to walk down the BST and add the new data as a node to the correct place"""

        if car.make == currentNode.make and car.model == currentNode.model:
            currentNode.cars.append(car)

        elif car < currentNode.cars[0]:                         # doesn't matter which car in the list we check
            if currentNode.getLeft():
                self._put(car, currentNode.getLeft())
            else:
                currentNode.setLeft(CarInventoryNode(car))
                currentNode.left.setParent(currentNode)

        else:
            if currentNode.getRight():
                self._put(car, currentNode.getRight())
            else:
                currentNode.setRight(CarInventoryNode(car))
                currentNode.right.setParent(currentNode)

    def doesCarExist(self, car):

        if self.root:
            node = self._get(car, self.root)

            if node:
                for resCar in node.cars:
                    if resCar == car:
                        return True

        return False

    def _get(self, car, currentNode):
        """Helper function to walk down the BST and get the node corresponding to the given make/model of given car"""

        if not currentNode:
            return None
        elif car.make == currentNode.make and car.model == currentNode.model:
            return currentNode
        elif car < currentNode.cars[0]:
            return self._get(car, currentNode.getLeft())
        else:
            return self._get(car, currentNode.getRight())

    def inOrder(self, currentNode="_"):

        if currentNode == "_":
            currentNode = self.root

        output = ""
        if currentNode:
            output += self.inOrder(currentNode.getLeft())
            output += str(currentNode)
            output += self.inOrder(currentNode.getRight())

        return output

    def preOrder(self, currentNode="_"):

        if currentNode == "_":              # Just a placeholder so I can use a recursive loop
            currentNode = self.root

        output = ""
        if currentNode:
            output += str(currentNode)
            output += self.preOrder(currentNode.getLeft())
            output += self.preOrder(currentNode.getRight())

        return output

    def postOrder(self, currentNode="_"):

        if currentNode == "_":
            currentNode = self.root

        output = ""
        if currentNode:
            output += self.postOrder(currentNode.getLeft())
            output += self.postOrder(currentNode.getRight())
            output += str(currentNode)

        return output

    def getBestCar(self, make, model):
        """Defined as the latest year and most expensive"""

        if self.root:
            tempCar = Car(make, model, 0, 0)                # Allows us to use previous helper function to find node
            node = self._get(tempCar, self.root)

            if node:

                best = node.cars[0]
                for resCar in node.cars:
                    if best < resCar:
                        best = resCar

                return best

        return None

    def getWorstCar(self, make, model):
        """Defined as the oldest year and the least expensive"""

        if self.root:
            tempCar = Car(make, model, 0, 0)
            node = self._get(tempCar, self.root)

            if node:

                worst = node.cars[0]
                for resCar in node.cars:
                    if worst > resCar:
                        worst = resCar

                return worst

        return None

    def getTotalInventoryPrice(self, currentNode="_"):

        if currentNode == "_":
            currentNode = self.root

        total = 0
        if currentNode:
            for car in currentNode.cars:
                total += car.price
            total += self.getTotalInventoryPrice(currentNode.getLeft())
            total += self.getTotalInventoryPrice(currentNode.getRight())

        return total

    def getSuccessor(self, make, model):

        if self.root:
            tempCar = Car(make, model, 0, 0)
            currentNode = self._get(tempCar, self.root)

            if currentNode:

                # Leftmost node in the right tree
                if currentNode.getRight():

                    successorNode = currentNode.getRight()
                    while successorNode.getLeft():
                        successorNode = successorNode.getLeft()

                    return successorNode

                # Parent with greater value than self (may not exist!)
                if currentNode.getParent():                         # Not useful for deleting but required nonetheless!

                    successorNode = currentNode.getParent()
                    while successorNode:
                        if successorNode.cars[0] > currentNode.cars[0]:       # Doesn't matter which car we check
                            return successorNode
                        successorNode = successorNode.getParent()

        return None

    def removeCar(self, make, model, year, price):

        if self.root:
            tempCar = Car(make, model, year, price)
            currentNode = self._get(tempCar, self.root)

            if currentNode:

                for index, car in enumerate(currentNode.cars):
                    if tempCar == car:
                        currentNode.cars.pop(index)
                        break

                if len(currentNode.cars) == 0:
                    self._remove(currentNode)

                return True

        return False

    def _remove(self, currentNode):

        # Case 1 : No children
        if not currentNode.getLeft() and not currentNode.getRight():
            if not currentNode.getParent():
                self.root = None
            elif currentNode == currentNode.getParent().getLeft():
                currentNode.getParent().setLeft(None)
            else:
                currentNode.getParent().setRight(None)

        # Case 3 : Node to remove has both children
        elif currentNode.getLeft() and currentNode.getRight():
            successorNode = self.getSuccessor(currentNode.make, currentNode.model)

            # Removing successor if no children
            if not successorNode.getLeft() and not successorNode.getRight():
                if successorNode == successorNode.getParent().getLeft():
                    successorNode.getParent().setLeft(None)
                else:
                    successorNode.getParent().setRight(None)

            # Removing successor if only right child (not possible to have left by definition)
            elif successorNode.getLeft() or successorNode.getRight():
                if successorNode.getRight():
                    if successorNode == successorNode.getParent().getLeft():
                        successorNode.getParent().setLeft(successorNode.getRight())
                    else:
                        successorNode.getParent().setRight(successorNode.getRight())
                    successorNode.getRight().setParent(successorNode.getParent())

            # Copying data from successor node to current node (not parent & children data as position is BST unchanged)
            currentNode.make = successorNode.make
            currentNode.model = successorNode.model
            currentNode.cars = successorNode.cars

        # Case 2 : Node to remove has one child
        else:

            # If the node to remove HAS a leftChild
            if currentNode.getLeft():
                if not currentNode.getParent():
                    # If current node is the root (no parent to reference to so set equal to its child)
                    currentNode.make = currentNode.getLeft().make
                    currentNode.model = currentNode.getLeft().model
                    currentNode.cars = currentNode.getLeft().cars
                    currentNode.setRight(currentNode.getLeft().getRight())      # Ordering of statements matters!!
                    # Not a problem in the lectures since the argument has already been passed through into function!
                    currentNode.setLeft(currentNode.getLeft().getLeft())
                    if currentNode.getLeft():
                        currentNode.getLeft().setParent(currentNode)
                    if currentNode.getRight():
                        currentNode.getRight().setParent(currentNode)
                elif currentNode == currentNode.getParent().getLeft():
                    currentNode.getLeft().setParent(currentNode.getParent())
                    currentNode.getParent().setLeft(currentNode.getLeft())
                elif currentNode == currentNode.getParent().getRight():
                    currentNode.getLeft().setParent(currentNode.getParent())
                    currentNode.getParent().setRight(currentNode.getLeft())

            # Similarly, if node to remove HAS a rightChild instead
            else:
                if not currentNode.getParent():
                    # If current node is the root (no parent to reference to so set equal to its child)
                    currentNode.make = currentNode.getRight().make
                    currentNode.model = currentNode.getRight().model
                    currentNode.cars = currentNode.getRight().cars
                    currentNode.setLeft(currentNode.getRight().getLeft())
                    currentNode.setRight(currentNode.getRight().getRight())
                    if currentNode.getLeft():
                        currentNode.getLeft().setParent(currentNode)
                    if currentNode.getRight():
                        currentNode.getRight().setParent(currentNode)
                elif currentNode == currentNode.getParent().getLeft():
                    currentNode.getRight().setParent(currentNode.getParent())
                    currentNode.getParent().setLeft(currentNode.getRight())
                elif currentNode == currentNode.getParent().getRight():
                    currentNode.getRight().setParent(currentNode.getParent())
                    currentNode.getParent().setRight(currentNode.getRight())


# In hindsight, better to work with .left, .right, .parent stuff as it's infinitely more clear and concise!
