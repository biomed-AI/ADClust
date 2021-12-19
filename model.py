import torch

class _ADClust_Autoencoder(torch.nn.Module):

    def __init__(self, input_dim: int, embedding_size: int):
        super(_ADClust_Autoencoder, self).__init__()

        self.encoder = torch.nn.Sequential(
            torch.nn.Linear(input_dim, 512),
            torch.nn.ReLU(inplace=True),
            torch.nn.Linear(512, 256),
            torch.nn.ReLU(inplace=True),
            torch.nn.Linear(256, 128),
            torch.nn.ReLU(inplace=True),
            torch.nn.Linear(128, embedding_size))

        self.decoder = torch.nn.Sequential(
            torch.nn.Linear(embedding_size, 128),
            torch.nn.ReLU(inplace=True),
            torch.nn.Linear(128, 256),
            torch.nn.ReLU(inplace=True),
            torch.nn.Linear(256, 512),
            torch.nn.ReLU(inplace=True),
            torch.nn.Linear(512, input_dim),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def decode(self, embedded: torch.Tensor) -> torch.Tensor:
        return self.decoder(embedded)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        embedded = self.encode(x)
        reconstruction = self.decode(embedded)
        return reconstruction

    def start_training(self, trainloader, n_epochs, device, optimizer, loss_fn):
        for _ in range(n_epochs):
            for batch, _ in trainloader:
                batch_data = batch.to(device)
                reconstruction = self.forward(batch_data)
                loss = loss_fn(reconstruction, batch_data)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
