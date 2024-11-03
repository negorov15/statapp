# Represents behaviour and logic.
# Needs sep. folder
class Model:

    def __init__(self, table, metadata):
        self.table = table
        self.metadata = metadata

    def access_database(self):
        ...

    def update_info(self, info):
        # updates info and sends back to controller
        ...

    def create(self, data):
        self.model.update_info(data)
        ...

    def delete(self, data):
        self.model.update_info(data)
        ...

    def edit(self, data):
        self.model.update_info(data)
        ...

    def transpose(self):
        transposed = ...
        self.model.update_info(transposed)
        ...