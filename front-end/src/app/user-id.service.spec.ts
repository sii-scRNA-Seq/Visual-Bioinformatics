import { HttpClientTestingModule } from '@angular/common/http/testing';
import { TestBed } from '@angular/core/testing';

import { UserIdService } from './user-id.service';
import { first, firstValueFrom } from 'rxjs';
import { BackendHttpClient } from './backend-http.client';
import { MatSnackBarModule } from '@angular/material/snack-bar';

describe('UserIdService', () => {
  let service: UserIdService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [ 
        HttpClientTestingModule,
        MatSnackBarModule,
       ]
    });
    service = TestBed.inject(UserIdService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('setUserId', () => {
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
    
    it('should update userId in local storage and in service as expected', async () => {
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
