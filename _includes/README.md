# Files in this folder and what they do:


- ## `site-header.html`
    `site-header.html` includes:

    - The `<head>` document information for the site

    - The sitewide header (including site title as determined in [`_config.yml`](../_config.yml))

    - The main navigation of the site, which will display only if the `displayMainNavigation` boolean in [`_data/siteNavigation.yml`](../_data/siteNavigation.yml) is set to `true`

    - Opening tags for `</main>` and `</body>`

- ## `site-footer.html`

    `site-footer.html` includes:

    - The sitewide footer and copyright information

    - Closing tags for `</main>` and `</body>`

    - Any scripts that are loaded before the closing `</body>` tag