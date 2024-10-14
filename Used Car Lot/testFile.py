from Car import Car
from CarInventoryNode import CarInventoryNode
from CarInventory import CarInventory

class TestCar:

    def test_constructor(self):

        car = Car("Dodge", "dart", 2015, 6000)
        assert car.make == 'DODGE'
        assert car.model == "DART"
        assert car.year == 2015
        assert car.price == 6000

    def test_inequalities(self):

        car1 = Car("Dodge", "dart", 2015, 6000)
        car2 = Car("Honda", "CRV", 2007, 8000)
        assert car2 > car1
        assert car1 < car2

        car2.make = 'dodge'.upper()                         # Not passing through the constructor so manually capitalize
        assert car1 > car2
        assert car2 < car1

        car2.model = 'Dart'.upper()
        assert car1 > car2
        assert car2 < car1

        car2.year = 2015
        assert car2 > car1
        assert car1 < car2

        car2.price = 6000
        assert car1 == car2

    def test_string(self):

        car1 = Car("Dodge", "dart", 2015, 6000)
        assert str(car1) == "Make: DODGE, Model: DART, Year: 2015, Price: $6000"

        car2 = Car("Honda", "CRV", 2007, 8000)
        assert str(car2) == "Make: HONDA, Model: CRV, Year: 2007, Price: $8000"


class TestCarInventoryNode:

    def test_addCar(self):

        car1 = Car("Dodge", "dart", 2015, 6000)
        car2 = Car("dodge", "DaRt", 2003, 5000)
        car3 = Car("dodge", "darT", 2008, 8000)

        carNode = CarInventoryNode(car1)
        carNode.cars.append(car2)
        assert str(carNode) == "Make: DODGE, Model: DART, Year: 2015, Price: $6000\n" \
                               "Make: DODGE, Model: DART, Year: 2003, Price: $5000\n"

        carNode.cars.append(car3)
        assert str(carNode) == "Make: DODGE, Model: DART, Year: 2015, Price: $6000\n" \
                               "Make: DODGE, Model: DART, Year: 2003, Price: $5000\n" \
                               "Make: DODGE, Model: DART, Year: 2008, Price: $8000\n"

        carNode = CarInventoryNode(car1)
        assert str(carNode) == str(car1) + "\n"


class TestCarInventory:

    def test_doesCarExist(self):

        bst = CarInventory()

        car1 = Car("Nissan", "Leaf", 2018, 18000)
        car2 = Car("Tesla", "Model3", 2018, 50000)
        car3 = Car("Dodge", "dart", 2015, 6000)
        car4 = Car("dodge", "DaRt", 2003, 5000)

        bst.addCar(car1)
        bst.addCar(car2)

        assert bst.doesCarExist(car2) is True
        assert bst.doesCarExist(car3) is False

        bst.addCar(car3)

        assert bst.doesCarExist(car3) is True
        assert bst.doesCarExist(car4) is False

        bst.addCar(car4)

        assert bst.doesCarExist(car4) is True

    def test_ordering(self):

        bst = CarInventory()
        car1 = Car("Nissan", "Leaf", 2018, 18000)
        car2 = Car("Tesla", "Model3", 2018, 50000)
        car3 = Car("Mercedes", "Sprinter", 2022, 40000)
        car4 = Car("Mercedes", "Sprinter", 2014, 25000)
        car5 = Car("Ford", "Ranger", 2021, 25000)
        bst.addCar(car1)
        bst.addCar(car2)
        bst.addCar(car3)
        bst.addCar(car4)
        bst.addCar(car5)

        assert bst.inOrder() == "Make: FORD, Model: RANGER, Year: 2021, Price: $25000\n" \
                                "Make: MERCEDES, Model: SPRINTER, Year: 2022, Price: $40000\n" \
                                "Make: MERCEDES, Model: SPRINTER, Year: 2014, Price: $25000\n" \
                                "Make: NISSAN, Model: LEAF, Year: 2018, Price: $18000\n" \
                                "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n"

        assert bst.preOrder() == "Make: NISSAN, Model: LEAF, Year: 2018, Price: $18000\n" \
                                 "Make: MERCEDES, Model: SPRINTER, Year: 2022, Price: $40000\n" \
                                 "Make: MERCEDES, Model: SPRINTER, Year: 2014, Price: $25000\n" \
                                 "Make: FORD, Model: RANGER, Year: 2021, Price: $25000\n" \
                                 "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n"

        assert bst.postOrder() == "Make: FORD, Model: RANGER, Year: 2021, Price: $25000\n" \
                                  "Make: MERCEDES, Model: SPRINTER, Year: 2022, Price: $40000\n" \
                                  "Make: MERCEDES, Model: SPRINTER, Year: 2014, Price: $25000\n" \
                                  "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n" \
                                  "Make: NISSAN, Model: LEAF, Year: 2018, Price: $18000\n"

        bst = CarInventory()
        assert bst.inOrder() == ""
        assert bst.preOrder() == ""
        assert bst.postOrder() == ""
        assert bst.doesCarExist(car4) is False

    def test_getBest_getWorst(self):

        bst = CarInventory()

        car1 = Car("Nissan", "Leaf", 2018, 18000)
        car2 = Car("Tesla", "Model3", 2018, 50000)
        car3 = Car("Mercedes", "Sprinter", 2022, 40000)
        car4 = Car("Mercedes", "Sprinter", 2014, 25000)
        car5 = Car("Ford", "Ranger", 2021, 25000)

        bst.addCar(car1)
        bst.addCar(car2)
        bst.addCar(car3)
        bst.addCar(car4)
        bst.addCar(car5)

        assert bst.getBestCar("Nissan", "Leaf") == car1
        assert bst.getBestCar("Mercedes", "Sprinter") == car3
        assert bst.getBestCar("Honda", "Accord") is None

        assert bst.getWorstCar("Nissan", "Leaf") == car1
        assert bst.getWorstCar("Mercedes", "Sprinter") == car4
        assert bst.getBestCar("Honda", "Accord") is None

        bst = CarInventory()
        assert bst.getBestCar("Nissan", "Leaf") is None
        assert bst.getWorstCar("Nissan", "Leaf") is None

    def test_getTotalInventoryPrice(self):

        bst = CarInventory()

        car1 = Car("Nissan", "Leaf", 2018, 18000)
        car2 = Car("Tesla", "Model3", 2018, 50000)
        car3 = Car("Mercedes", "Sprinter", 2022, 40000)
        car4 = Car("Mercedes", "Sprinter", 2014, 25000)
        car5 = Car("Ford", "Ranger", 2021, 25000)

        bst.addCar(car1)
        bst.addCar(car2)
        bst.addCar(car3)
        bst.addCar(car4)
        bst.addCar(car5)

        assert bst.getTotalInventoryPrice() == 158000

        bst = CarInventory()
        assert bst.getTotalInventoryPrice() == 0

    def test_getSuccessor(self):

        bst = CarInventory()

        car1 = Car("Nissan", "Leaf", 2018, 18000)
        car2 = Car("Tesla", "Model3", 2018, 50000)
        car3 = Car("Mercedes", "Sprinter", 2022, 40000)
        car4 = Car("Mercedes", "Sprinter", 2014, 25000)
        car5 = Car("Ford", "Ranger", 2021, 25000)
        car6 = Car("Mazda", "CX-5", 2022, 25000)
        car7 = Car("Tesla", "Model3", 2018, 50000)
        car8 = Car("BMW", "X5", 2022, 60000)
        car9 = Car("BMW", "X5", 2020, 58000)
        car10 = Car("Audi", "A3", 2021, 25000)

        bst.addCar(car1)
        bst.addCar(car2)
        bst.addCar(car3)
        bst.addCar(car4)
        bst.addCar(car5)
        bst.addCar(car6)
        bst.addCar(car7)
        bst.addCar(car8)
        bst.addCar(car9)
        bst.addCar(car10)

        assert str(bst.getSuccessor("Ford", "Ranger")) == "Make: MAZDA, Model: CX-5, Year: 2022, Price: $25000\n"
        assert str(bst.getSuccessor("BMW", "X5")) == "Make: FORD, Model: RANGER, Year: 2021, Price: $25000\n"
        assert bst.getSuccessor("Tesla", "Model3") is None
        assert str(bst.getSuccessor("Mazda", "CX-5")) == "Make: MERCEDES, Model: SPRINTER, Year: 2022, Price: $40000\n" \
                                                         "Make: MERCEDES, Model: SPRINTER, Year: 2014, Price: $25000\n"

    def test_getSuccessor_edgeCases(self):

        bst = CarInventory()

        car1 = Car("Nissan", "Leaf", 2018, 18000)
        car3 = Car("Mercedes", "Sprinter", 2022, 40000)
        car5 = Car("Ford", "Ranger", 2021, 25000)
        car6 = Car("Mazda", "CX-5", 2022, 25000)
        car9 = Car("BMW", "X5", 2020, 58000)
        car10 = Car("Audi", "A3", 2021, 25000)

        bst.addCar(car1)
        bst.addCar(car5)
        bst.addCar(car9)
        bst.addCar(car10)
        bst.addCar(car6)
        bst.addCar(car3)

        assert bst.getSuccessor("Nissan", "Leaf") is None
        assert str(bst.getSuccessor("BMW", "X5")) == "Make: FORD, Model: RANGER, Year: 2021, Price: $25000\n"
        assert str(bst.getSuccessor("Mazda", "CX-5")) == "Make: MERCEDES, Model: SPRINTER, Year: 2022, Price: $40000\n"
        assert str(bst.getSuccessor("Mercedes", "Sprinter")) == "Make: NISSAN, Model: LEAF, Year: 2018, Price: $18000\n"

    def test_removeCar(self):

        bst = CarInventory()

        car1 = Car("Mazda", "CX-5", 2022, 25000)
        car2 = Car("Tesla", "Model3", 2018, 50000)
        car3 = Car("BMW", "X5", 2022, 60000)
        car4 = Car("BMW", "X5", 2020, 58000)
        car5 = Car("Audi", "A3", 2021, 25000)

        bst.addCar(car1)
        bst.addCar(car2)
        bst.addCar(car3)
        bst.addCar(car4)
        bst.addCar(car5)

        bst.removeCar("BMW", "X5", 2020, 58000)
        assert bst.inOrder() == "Make: AUDI, Model: A3, Year: 2021, Price: $25000\n" \
                                "Make: BMW, Model: X5, Year: 2022, Price: $60000\n" \
                                "Make: MAZDA, Model: CX-5, Year: 2022, Price: $25000\n" \
                                "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n"
        assert bst.preOrder() == "Make: MAZDA, Model: CX-5, Year: 2022, Price: $25000\n" \
                                 "Make: BMW, Model: X5, Year: 2022, Price: $60000\n" \
                                 "Make: AUDI, Model: A3, Year: 2021, Price: $25000\n" \
                                 "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n"
        assert bst.postOrder() == "Make: AUDI, Model: A3, Year: 2021, Price: $25000\n" \
                                  "Make: BMW, Model: X5, Year: 2022, Price: $60000\n" \
                                  "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n" \
                                  "Make: MAZDA, Model: CX-5, Year: 2022, Price: $25000\n"

        bst.removeCar("BMW", "X5", 2022, 60000)
        assert bst.inOrder() == "Make: AUDI, Model: A3, Year: 2021, Price: $25000\n" \
                                "Make: MAZDA, Model: CX-5, Year: 2022, Price: $25000\n" \
                                "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n"
        assert bst.preOrder() == "Make: MAZDA, Model: CX-5, Year: 2022, Price: $25000\n" \
                                 "Make: AUDI, Model: A3, Year: 2021, Price: $25000\n" \
                                 "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n"
        assert bst.postOrder() == "Make: AUDI, Model: A3, Year: 2021, Price: $25000\n" \
                                  "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n" \
                                  "Make: MAZDA, Model: CX-5, Year: 2022, Price: $25000\n"

    def test_removeCar_edgeCases(self):

        bst = CarInventory()

        car1 = Car("Mazda", "CX-5", 2022, 25000)
        car2 = Car("Tesla", "Model3", 2018, 50000)
        car3 = Car("BMW", "X5", 2022, 60000)

        bst.addCar(car1)
        bst.addCar(car3)

        bst.removeCar("Mazda", "CX-5", 2022, 25000)
        assert bst.preOrder() == "Make: BMW, Model: X5, Year: 2022, Price: $60000\n"

        bst = CarInventory()
        bst.addCar(car1)
        bst.addCar(car2)

        bst.removeCar("Mazda", "CX-5", 2022, 25000)
        assert bst.preOrder() == "Make: TESLA, Model: MODEL3, Year: 2018, Price: $50000\n"

        bst.removeCar("Tesla", "Model3", 2018, 50000)
        assert bst.inOrder() == ""


