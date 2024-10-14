class Car:

    def __init__(self, make: str, model: str, year: int, price: int):
        self.make = make.upper()
        self.model = model.upper()
        self.year = year
        self.price = price

    def __gt__(self, other):

        if self.make == other.make:
            if self.model == other.model:
                if self.year == other.year:
                    return self.price > other.price
                else:
                    return self.year > other.year
            else:
                return self.model > other.model
        else:
            return self.make > other.make

        # Returns false by default for an equality of two equal car objects

    def __lt__(self, other):

        if self.make == other.make:
            if self.model == other.model:
                if self.year == other.year:
                    return self.price < other.price
                else:
                    return self.year < other.year
            else:
                return self.model < other.model
        else:
            return self.make < other.make

        # Only reason this need to be 'duplicated' is so that for an equality it returns false like above not True?

    def __eq__(self, other):

        if self.make == other.make and self.model == other.model and self.year == other.year and \
                self.price == other.price:
            return True

        return False

    def __str__(self):
        return f"Make: {self.make}, Model: {self.model}, Year: {self.year}, Price: ${self.price}"
