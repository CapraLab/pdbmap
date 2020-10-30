# Make sure all the tests the config file that provides SQL and filesystem access to resources
import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--config", action="store", 
        default="../../config/global.config", 
        help="A -c option is required to load a db and filesystem configuration"
        )

@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--config")
