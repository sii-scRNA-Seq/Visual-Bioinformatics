import { MatSnackBarModule } from '@angular/material/snack-bar';
import { TestBed } from '@angular/core/testing';
import { first, firstValueFrom } from 'rxjs';
import { provideHttpClient, withInterceptorsFromDi } from '@angular/common/http';
import { provideHttpClientTesting } from '@angular/common/http/testing';

import { BackendHttpClient } from './backend-http.client';
import { UserIdService } from './user-id.service';

describe('UserIdService', () => {
  let service: UserIdService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        MatSnackBarModule
      ],
      providers: [
        provideHttpClient(withInterceptorsFromDi()), provideHttpClientTesting()
      ]
    });
  });

  it('should be created', () => {
    service = TestBed.inject(UserIdService);
    expect(service).toBeTruthy();
  });

  describe('constructor', () => {
    it('should result in datasetInfo being loaded', async () => {
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      spyOn(backendHttpClient, 'getUserId').and.returnValue(Promise.resolve('mock_user_id'));
      service = await TestBed.inject(UserIdService);
      const userId = await firstValueFrom(service.userId.pipe(first()));
      expect(userId).toEqual('mock_user_id');
    });
  });

  describe('setUserId', () => {
    beforeEach(() => {
      service = TestBed.inject(UserIdService);
    });

    it('should call client with local storage userId if it exists', async () => {
      spyOn(localStorage, 'getItem').and.returnValue('mock_user_id');
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      const backendHttpClientSpy = spyOn(backendHttpClient, 'getUserId');
      await service.setUserId();
      expect(backendHttpClientSpy).toHaveBeenCalledOnceWith('mock_user_id');
    });

    it('should call client with empty string if local storage userId does not exist', async () => {
      spyOn(localStorage, 'getItem').and.returnValue(null);
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      const backendHttpClientSpy = spyOn(backendHttpClient, 'getUserId');
      await service.setUserId();
      expect(backendHttpClientSpy).toHaveBeenCalledOnceWith('');
    });

    it('should update userId in local storage and in attribute as expected', async () => {
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      spyOn(backendHttpClient, 'getUserId').and.returnValue(Promise.resolve('mock_user_id'));
      const localStorageSpy = spyOn(localStorage, 'setItem');
      await service.setUserId();
      expect(localStorageSpy).toHaveBeenCalledOnceWith('userId', 'mock_user_id');
      const response = await firstValueFrom(service.userId.pipe(first()));
      expect(response).toEqual('mock_user_id');
    });
  });
});
