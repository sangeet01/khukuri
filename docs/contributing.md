# Contributing Guide

## Development Setup

```bash
git clone <repository>
cd khukuri
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

## Code Standards

### Style Guide
- Follow PEP 8
- Use Black for formatting: `black src/`
- Use flake8 for linting: `flake8 src/`
- Type hints required for public APIs

### Documentation
- Docstrings for all public functions/classes
- Google-style docstring format
- Update API reference when adding features

### Testing
- Write tests for all new features
- Maintain >80% code coverage
- Run tests before committing: `pytest tests/`

## Workflow

1. **Fork** the repository
2. **Create branch**: `git checkout -b feature/my-feature`
3. **Make changes** following code standards
4. **Write tests** for new functionality
5. **Run tests**: `pytest tests/ -v`
6. **Format code**: `black src/`
7. **Commit**: `git commit -m "Add feature X"`
8. **Push**: `git push origin feature/my-feature`
9. **Create Pull Request**

## Adding New Modules

1. Create directory: `src/my_module/`
2. Add `__init__.py` with exports
3. Implement functionality in separate files
4. Add tests in `tests/test_my_module/`
5. Update documentation
6. Add to `REPO_STRUCTURE.md`

## Running Tests

```bash
# All tests
pytest tests/ -v

# Specific module
pytest tests/test_core/ -v

# With coverage
pytest tests/ --cov=src --cov-report=html

# Skip slow tests
pytest tests/ -m "not slow"
```

## Code Review Checklist

- [ ] Code follows PEP 8
- [ ] All tests pass
- [ ] New tests added for new features
- [ ] Documentation updated
- [ ] No mock implementations
- [ ] Error handling implemented
- [ ] Logging added
- [ ] Type hints included

## Questions?

Open an issue for discussion before major changes.
