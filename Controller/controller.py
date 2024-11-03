class Controller:

    def __init__(self, model, view):
        self.model = model
        self.view = view

    # Gives some tasks to the model. For example: insert, create, delete
    def manipulate(self, user_input):
        self.model.create(user_input)
        ...

    # Passes the data to the view to generate an image.
    def show(self, data):
        self.view.show(data)
        ...