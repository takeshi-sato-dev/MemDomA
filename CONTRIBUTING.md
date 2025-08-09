# Contributing to MemDomA

We welcome contributions to MemDomA! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/MemDomA.git
   cd MemDomA/memdoma
   ```
3. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature-name
   ```

## Development Setup

1. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

2. Install the package in development mode:
   ```bash
   pip install -r requirements.txt
   pip install -e .
   ```

3. Install development dependencies:
   ```bash
   pip install pytest pytest-cov black flake8
   ```

## Testing Your Changes

Before submitting a pull request, please run the following tests:

### 1. Basic Import Test
```bash
python -c "from memdoma_package import MemDomA; print('✓ Import successful')"
```

### 2. Run Test with Sample Data
```bash
python examples/test_with_sample_data.py
```

### 3. Unit Tests (if available)
```bash
pytest tests/
```

### 4. Check Code Style
```bash
# Format code with black
black memdoma_package/

# Check with flake8
flake8 memdoma_package/ --max-line-length=100
```

### 5. Full Test Before PR

Run this checklist before submitting:

- [ ] Code runs without errors
- [ ] Test script passes (`python examples/test_with_sample_data.py`)
- [ ] Code follows PEP 8 style guidelines
- [ ] Documentation/docstrings updated if needed
- [ ] New features have example usage
- [ ] Commit messages are clear and descriptive

## Code Style

- Follow PEP 8
- Use descriptive variable names
- Add docstrings to all functions and classes
- Keep functions focused and small
- Comment complex logic

## Submitting Changes

1. Commit your changes:
   ```bash
   git add .
   git commit -m "Clear description of changes"
   ```

2. Push to your fork:
   ```bash
   git push origin feature-name
   ```

3. Submit a pull request through GitHub

## Pull Request Guidelines

- Describe what changes you made and why
- Reference any related issues
- Include example usage if adding new features
- Ensure all tests pass
- Update documentation if necessary

## Reporting Issues

When reporting issues, please include:
- Your operating system and Python version
- MemDomA version
- Complete error messages
- Minimal code example to reproduce the issue
- Input file specifications (number of atoms, frames, etc.)

## Feature Requests

We welcome feature requests! Please:
- Check existing issues first
- Describe the feature clearly
- Explain use cases
- Consider submitting a PR if you can implement it

## Questions

For questions about using MemDomA, please:
- Check the README and documentation first
- Look through existing GitHub issues
- Open a new issue with the "question" label

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to MemDomA!