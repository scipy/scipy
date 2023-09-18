# SPDX-License-Identifier: MIT

from __future__ import annotations

import collections
import dataclasses
import os
import os.path
import pathlib
import re
import typing

from collections.abc import Mapping
from typing import Any

import packaging.markers
import packaging.requirements
import packaging.specifiers
import packaging.version


__version__ = '0.7.1'


class ConfigurationError(Exception):
    '''Error in the backend metadata.'''
    def __init__(self, msg: str, *, key: str | None = None):
        super().__init__(msg)
        self._key = key

    @property
    def key(self) -> str | None:  # pragma: no cover
        return self._key


class RFC822Message():
    '''Python-flavored RFC 822 message implementation.'''

    def __init__(self) -> None:
        self.headers: collections.OrderedDict[str, list[str]] = collections.OrderedDict()
        self.body: str | None = None

    def __setitem__(self, name: str, value: str | None) -> None:
        if not value:
            return
        if name not in self.headers:
            self.headers[name] = []
        self.headers[name].append(value)

    def __str__(self) -> str:
        text = ''
        for name, entries in self.headers.items():
            for entry in entries:
                lines = entry.strip('\n').split('\n')
                text += f'{name}: {lines[0]}\n'
                for line in lines[1:]:
                    text += ' ' * 8 + line + '\n'
        if self.body:
            text += '\n' + self.body
        return text

    def __bytes__(self) -> bytes:
        return str(self).encode()


class DataFetcher():
    def __init__(self, data: Mapping[str, Any]) -> None:
        self._data = data

    def __contains__(self, key: Any) -> bool:
        if not isinstance(key, str):
            return False
        val = self._data
        try:
            for part in key.split('.'):
                val = val[part]
        except KeyError:
            return False
        return True

    def get(self, key: str) -> Any:
        val = self._data
        for part in key.split('.'):
            val = val[part]
        return val

    def get_str(self, key: str) -> str | None:
        try:
            val = self.get(key)
            if not isinstance(val, str):
                raise ConfigurationError(
                    f'Field `{key}` has an invalid type, '
                    f'expecting a string (got `{val}`)',
                    key=key,
                )
            return val
        except KeyError:
            return None

    def get_list(self, key: str) -> list[str]:
        try:
            val = self.get(key)
            if not isinstance(val, list):
                raise ConfigurationError(
                    f'Field `{key}` has an invalid type, '
                    f'expecting a list of strings (got `{val}`)',
                    key=val,
                )
            for item in val:
                if not isinstance(item, str):
                    raise ConfigurationError(
                        f'Field `{key}` contains item with invalid type, '
                        f'expecting a string (got `{item}`)',
                        key=key,
                    )
            return val
        except KeyError:
            return []

    def get_dict(self, key: str) -> dict[str, str]:
        try:
            val = self.get(key)
            if not isinstance(val, dict):
                raise ConfigurationError(
                    f'Field `{key}` has an invalid type, '
                    f'expecting a dictionary of strings (got `{val}`)',
                    key=key,
                )
            for subkey, item in val.items():
                if not isinstance(item, str):
                    raise ConfigurationError(
                        f'Field `{key}.{subkey}` has an invalid type, '
                        f'expecting a string (got `{item}`)',
                        key=f'{key}.{subkey}',
                    )
            return val
        except KeyError:
            return {}

    def get_people(self, key: str) -> list[tuple[str, str]]:
        try:
            val = self.get(key)
            if not (
                isinstance(val, list)
                and all(isinstance(x, dict) for x in val)
                and all(
                    isinstance(item, str)
                    for items in [_dict.values() for _dict in val]
                    for item in items
                )
            ):
                raise ConfigurationError(
                    f'Field `{key}` has an invalid type, expecting a list of '
                    f'dictionaries containing the `name` and/or `email` keys (got `{val}`)',
                    key=key,
                )
            return [
                (entry.get('name', 'Unknown'), entry.get('email'))
                for entry in val
            ]
        except KeyError:
            return []


class License(typing.NamedTuple):
    text: str
    file: pathlib.Path | None


class Readme(typing.NamedTuple):
    text: str
    file: pathlib.Path | None
    content_type: str


@dataclasses.dataclass
class StandardMetadata():
    name: str
    version: packaging.version.Version | None = None
    description: str | None = None
    license: License | None = None
    readme: Readme | None = None
    requires_python: packaging.specifiers.SpecifierSet | None = None
    dependencies: list[packaging.requirements.Requirement] = dataclasses.field(default_factory=list)
    optional_dependencies: dict[str, list[packaging.requirements.Requirement]] = dataclasses.field(default_factory=dict)
    entrypoints: dict[str, dict[str, str]] = dataclasses.field(default_factory=dict)
    authors: list[tuple[str, str]] = dataclasses.field(default_factory=list)
    maintainers: list[tuple[str, str]] = dataclasses.field(default_factory=list)
    urls: dict[str, str] = dataclasses.field(default_factory=dict)
    classifiers: list[str] = dataclasses.field(default_factory=list)
    keywords: list[str] = dataclasses.field(default_factory=list)
    scripts: dict[str, str] = dataclasses.field(default_factory=dict)
    gui_scripts: dict[str, str] = dataclasses.field(default_factory=dict)
    dynamic: list[str] = dataclasses.field(default_factory=list)

    def __post_init__(self) -> None:
        self.name = re.sub(r'[-_.]+', '-', self.name).lower()
        self._update_dynamic(self.version)

    @classmethod
    def from_pyproject(
        cls,
        data: Mapping[str, Any],
        project_dir: str | os.PathLike[str] = os.path.curdir,
    ) -> StandardMetadata:
        fetcher = DataFetcher(data)
        project_dir = pathlib.Path(project_dir)

        if 'project' not in fetcher:
            raise ConfigurationError('Section `project` missing in pyproject.toml')

        dynamic = fetcher.get_list('project.dynamic')
        if 'name' in dynamic:
            raise ConfigurationError('Unsupported field `name` in `project.dynamic`')

        for field in dynamic:
            if field in data['project']:
                raise ConfigurationError(
                    f'Field `project.{field}` declared as dynamic in but is defined'
                )

        name = fetcher.get_str('project.name')
        if not name:
            raise ConfigurationError('Field `project.name` missing')

        version_string = fetcher.get_str('project.version')
        requires_python_string = fetcher.get_str('project.requires-python')

        return cls(
            name,
            packaging.version.Version(version_string) if version_string else None,
            fetcher.get_str('project.description'),
            cls._get_license(fetcher, project_dir),
            cls._get_readme(fetcher, project_dir),
            packaging.specifiers.SpecifierSet(requires_python_string) if requires_python_string else None,
            cls._get_dependencies(fetcher),
            cls._get_optional_dependencies(fetcher),
            cls._get_entrypoints(fetcher),
            fetcher.get_people('project.authors'),
            fetcher.get_people('project.maintainers'),
            fetcher.get_dict('project.urls'),
            fetcher.get_list('project.classifiers'),
            fetcher.get_list('project.keywords'),
            fetcher.get_dict('project.scripts'),
            fetcher.get_dict('project.gui-scripts'),
            dynamic,
        )

    def _update_dynamic(self, value: Any) -> None:
        if value and 'version' in self.dynamic:
            self.dynamic.remove('version')
        elif not value and 'version' not in self.dynamic:
            self.dynamic.append('version')

    def __setattr__(self, name: str, value: Any) -> None:
        # update dynamic when version is set
        if name == 'version' and hasattr(self, 'dynamic'):
            self._update_dynamic(value)
        super().__setattr__(name, value)

    def as_rfc822(self) -> RFC822Message:
        message = RFC822Message()
        self.write_to_rfc822(message)
        return message

    def write_to_rfc822(self, message: RFC822Message) -> None:  # noqa: C901
        message['Metadata-Version'] = '2.2' if self.dynamic else '2.1'
        message['Name'] = self.name
        if not self.version:
            raise ConfigurationError('Missing version field')
        message['Version'] = str(self.version)
        # skip 'Platform'
        # skip 'Supported-Platform'
        if self.description:
            message['Summary'] = self.description
        message['Keywords'] = ' '.join(self.keywords)
        if 'homepage' in self.urls:
            message['Home-page'] = self.urls['homepage']
        # skip 'Download-URL'
        message['Author'] = self._name_list(self.authors)
        message['Author-Email'] = self._email_list(self.authors)
        message['Maintainer'] = self._name_list(self.maintainers)
        message['Maintainer-Email'] = self._email_list(self.maintainers)
        if self.license:
            message['License'] = self.license.text
        for classifier in self.classifiers:
            message['Classifier'] = classifier
        # skip 'Provides-Dist'
        # skip 'Obsoletes-Dist'
        # skip 'Requires-External'
        for name, url in self.urls.items():
            message['Project-URL'] = f'{name.capitalize()}, {url}'
        if self.requires_python:
            message['Requires-Python'] = str(self.requires_python)
        for dep in self.dependencies:
            message['Requires-Dist'] = str(dep)
        for extra, requirements in self.optional_dependencies.items():
            message['Provides-Extra'] = extra
            for requirement in requirements:
                message['Requires-Dist'] = str(self._build_extra_req(extra, requirement))
        if self.readme:
            if self.readme.content_type:
                message['Description-Content-Type'] = self.readme.content_type
            message.body = self.readme.text
        # Core Metadata 2.2
        for field in self.dynamic:
            if field in ('name', 'version'):
                raise ConfigurationError(f'Field cannot be dynamic: {field}')
            message['Dynamic'] = field

    def _name_list(self, people: list[tuple[str, str]]) -> str:
        return ', '.join(
            name
            for name, email_ in people
            if not email_
        )

    def _email_list(self, people: list[tuple[str, str]]) -> str:
        return ', '.join([
            '{}{}'.format(name, f' <{_email}>' if _email else '')
            for name, _email in people
            if _email
        ])

    def _build_extra_req(
        self,
        extra: str,
        requirement: packaging.requirements.Requirement,
    ) -> packaging.requirements.Requirement:
        if requirement.marker:  # append our extra to the marker
            requirement.marker = packaging.markers.Marker(
                str(requirement.marker) + f' and extra == "{extra}"'
            )
        else:  # add our extra marker
            requirement.marker = packaging.markers.Marker(f'extra == "{extra}"')
        return requirement

    @staticmethod
    def _get_license(fetcher: DataFetcher, project_dir: pathlib.Path) -> License | None:
        if 'project.license' not in fetcher:
            return None

        _license = fetcher.get_dict('project.license')
        for field in _license:
            if field not in ('file', 'text'):
                raise ConfigurationError(
                    f'Unexpected field `project.license.{field}`',
                    key=f'project.license.{field}',
                )

        file: pathlib.Path | None = None
        filename = fetcher.get_str('project.license.file')
        text = fetcher.get_str('project.license.text')

        if (filename and text) or (not filename and not text):
            raise ConfigurationError(
                f'Invalid `project.license` value, expecting either `file` or `text` (got `{_license}`)',
                key='project.license',
            )

        if filename:
            file = project_dir.joinpath(filename)
            if not file.is_file():
                raise ConfigurationError(
                    f'License file not found (`{filename}`)',
                    key='project.license.file',
                )
            text = file.read_text(encoding='utf-8')

        assert text is not None
        return License(text, file)

    @staticmethod
    def _get_readme(fetcher: DataFetcher, project_dir: pathlib.Path) -> Readme | None:  # noqa: C901
        if 'project.readme' not in fetcher:
            return None

        filename: str | None
        file: pathlib.Path | None = None
        text: str | None
        content_type: str | None

        readme = fetcher.get('project.readme')
        if isinstance(readme, str):
            # readme is a file
            text = None
            filename = readme
            if filename.endswith('.md'):
                content_type = 'text/markdown'
            elif filename.endswith('.rst'):
                content_type = 'text/x-rst'
            else:
                raise ConfigurationError(
                    f'Could not infer content type for readme file `{filename}`',
                    key='project.readme',
                )
        elif isinstance(readme, dict):
            # readme is a dict containing either 'file' or 'text', and content-type
            for field in readme:
                if field not in ('content-type', 'file', 'text'):
                    raise ConfigurationError(
                        f'Unexpected field `project.readme.{field}`',
                        key=f'project.readme.{field}',
                    )
            content_type = fetcher.get_str('project.readme.content-type')
            filename = fetcher.get_str('project.readme.file')
            text = fetcher.get_str('project.readme.text')
            if (filename and text) or (not filename and not text):
                raise ConfigurationError(
                    f'Invalid `project.readme` value, expecting either `file` or `text` (got `{readme}`)',
                    key='project.license',
                )
            if not content_type:
                raise ConfigurationError(
                    'Field `project.readme.content-type` missing',
                    key='project.readme.content-type',
                )
        else:
            raise ConfigurationError(
                f'Field `project.readme` has an invalid type, expecting either, '
                f'a string or dictionary of strings (got `{readme}`)',
                key='project.readme',
            )

        if filename:
            file = project_dir.joinpath(filename)
            if not file.is_file():
                raise ConfigurationError(
                    f'Readme file not found (`{filename}`)',
                    key='project.license.file',
                )
            text = file.read_text(encoding='utf-8')

        assert text is not None
        return Readme(text, file, content_type)

    @staticmethod
    def _get_dependencies(fetcher: DataFetcher) -> list[packaging.requirements.Requirement]:
        try:
            requirement_strings = fetcher.get_list('project.dependencies')
        except KeyError:
            return []

        requirements: list[packaging.requirements.Requirement] = []
        for req in requirement_strings:
            try:
                requirements.append(packaging.requirements.Requirement(req))
            except packaging.requirements.InvalidRequirement as e:
                raise ConfigurationError(
                    'Field `project.dependencies` contains an invalid PEP 508 '
                    f'requirement string `{req}` (`{str(e)}`)'
                )
        return requirements

    @staticmethod
    def _get_optional_dependencies(fetcher: DataFetcher) -> dict[str, list[packaging.requirements.Requirement]]:
        try:
            val = fetcher.get('project.optional-dependencies')
        except KeyError:
            return {}

        requirements_dict: collections.defaultdict[str, list[packaging.requirements.Requirement]] = collections.defaultdict(list)
        if not isinstance(val, dict):
            raise ConfigurationError(
                'Field `project.optional-dependencies` has an invalid type, expecting a '
                f'dictionary of PEP 508 requirement strings (got `{val}`)'
            )
        for extra, requirements in val.copy().items():
            assert isinstance(extra, str)
            if not isinstance(requirements, list):
                raise ConfigurationError(
                    f'Field `project.optional-dependencies.{extra}` has an invalid type, expecting a '
                    f'dictionary PEP 508 requirement strings (got `{requirements}`)'
                )
            for i, req in enumerate(requirements):
                if not isinstance(req, str):
                    raise ConfigurationError(
                        f'Field `project.optional-dependencies.{extra}` has an invalid type, '
                        f'expecting a PEP 508 requirement string (got `{req}`)'
                    )
                try:
                    requirements_dict[extra].append(packaging.requirements.Requirement(req))
                except packaging.requirements.InvalidRequirement as e:
                    raise ConfigurationError(
                        f'Field `project.optional-dependencies.{extra}` contains '
                        f'an invalid PEP 508 requirement string `{req}` (`{str(e)}`)'
                    )
        return dict(requirements_dict)

    @staticmethod
    def _get_entrypoints(fetcher: DataFetcher) -> dict[str, dict[str, str]]:
        try:
            val = fetcher.get('project.entry-points')
        except KeyError:
            return {}
        if not isinstance(val, dict):
            raise ConfigurationError(
                'Field `project.entry-points` has an invalid type, expecting a '
                f'dictionary of entrypoint sections (got `{val}`)'
            )
        for section, entrypoints in val.items():
            assert isinstance(section, str)
            if not isinstance(entrypoints, dict):
                raise ConfigurationError(
                    f'Field `project.entry-points.{section}` has an invalid type, expecting a '
                    f'dictionary of entrypoints (got `{entrypoints}`)'
                )
            for name, entrypoint in entrypoints.items():
                assert isinstance(name, str)
                if not isinstance(entrypoint, str):
                    raise ConfigurationError(
                        f'Field `project.entry-points.{section}.{name}` has an invalid type, '
                        f'expecting a string (got `{entrypoint}`)'
                    )
        return val
