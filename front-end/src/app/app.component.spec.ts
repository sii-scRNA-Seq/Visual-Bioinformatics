import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';
import { TestBed } from '@angular/core/testing';
import { provideHttpClient, withInterceptorsFromDi } from '@angular/common/http';
import { provideHttpClientTesting } from '@angular/common/http/testing';

import { AppComponent } from './app.component';
import { BlockLibraryComponent } from './block-library/block-library.component';
import { CanvasComponent } from './canvas/canvas.component';
import { MockSocket } from './mock-socket';
import { OutputDisplayComponent } from './output-display/output-display.component';
import { SOCKET } from './socket';

describe('AppComponent', () => {
  beforeEach(() => TestBed.configureTestingModule({
    declarations: [
      AppComponent,
      BlockLibraryComponent,
      CanvasComponent,
      OutputDisplayComponent
    ],
    imports: [
      MatCardModule,
      MatSnackBarModule
    ],
    providers: [
      { provide: SOCKET, useClass: MockSocket },
      provideHttpClient(withInterceptorsFromDi()),
      provideHttpClientTesting()
    ]
  }));

  it('should create the app', () => {
    const fixture = TestBed.createComponent(AppComponent);
    const app = fixture.componentInstance;
    expect(app).toBeTruthy();
  });
});
